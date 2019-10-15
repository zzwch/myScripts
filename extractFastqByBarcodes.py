from Bio.SeqIO.QualityIO import FastqGeneralIterator
import io, os, time, gzip
import multiprocessing
import logging
import click

def pairedFastqExtractor(fq1, fq2, fq, outdir, barcodes, mismatch, logging):
    #buffer_max = 100000000
    mismatch = int(mismatch) # in case of str
    barcodes_mis_dict = mismatch_dict(barcodes, mismatch)
    #out1 = io.BufferedWriter(gzip.open(fq, 'w'), buffer_size = buffer_max)
    out1 = gzip.open(os.path.join(outdir, os.path.basename(fq)+'_R1.fastq.gz'), 'wt')
    out2 = gzip.open(os.path.join(outdir, os.path.basename(fq)+'_R2.fastq.gz'), 'wt')
    #bbcount = dict(zip(barcodes, [[[0 for i in range(mismatch+1)] for j in range(2)] for k in range(len(barcodes))]))
    #bbcount['ambiguous'] = [0 for i in range(2)]
    #bbcount['unmatched'] = [0 for j in range(2)]
    logging.info('Start processing ' + fq + '...')
    if len(fq1) == len(fq2):
        for i in range(0, len(fq1)):
            #in1 = io.BufferedReader(gzip.open(fq1[i],'r'), buffer_size = buffer_max)
            #in2 = io.BufferedReader(gzip.open(fq2[i],'r'), buffer_size = buffer_max)
            in1 = gzip.open(fq1[i],'rt')
            in2 = gzip.open(fq2[i],'rt')
            r1 = FastqGeneralIterator(in1)
            r2 = FastqGeneralIterator(in2)
            NR = 0
            count = 0
            for info in r2:
                seq = r1.__next__()
                tag = info[1][0:8]
                #isbar_flag = False
                #if tag in barcodes_mis_dict:
                    #bb = barcodes_mis_dict[tag]
                    #bbcount[bb[1]][0][bb[0]] += 1
                #    isbar_flag = True
                #else:
                #    bbcount['unmatched'][0] += 1
                #if isbar_flag:
                if tag in barcodes_mis_dict:
                    count += 1
                    out1.write('@%s\n%s\n+\n%s\n' % seq)
                    out2.write('@%s\n%s\n+\n%s\n' % info)
                    #bbcount[bb[1]][1][bb[0]] += 1
                #else:
                    #bbcount['unmatched'][1] += 1
                NR += 1
                if NR % 100000 == 0:
                    logging.info('{} (run_{}): processed {}/{}'.format(fq, str(i+1), str(count), str(NR)))
            logging.info('{} (run_{}): Done {}/{}. \n'.format(fq, str(i+1), str(count), str(NR)))
            in1.close()
            in2.close()
    else:
        return None
    out1.flush()
    out2.flush()
    out1.close()
    out2.close()
    return None

def mismatch_dict(barcodes, mismatch):
    # using hash to save time, Great time when iterations heavily used
    barcodes_mis_dict = {}
    for bar in barcodes:
        if mismatch == 0:
            barcodes_mis_dict[bar] = (0, bar)
        elif mismatch == 1:
            for i in range(len(bar)):
                for b in 'ATGCN':
                    bar_mis = bar[0:i] + b + bar[i+1:]
                    barcodes_mis_dict[bar_mis] = (hamming2(bar_mis, bar), bar)
        elif mismatch == 2:
            for i in range(len(bar)-1):
                for j in range(i+1, len(bar)):
                    for b1 in 'ATGCN':
                        for b2 in 'ATGCN':
                            bar_mis = bar[0:i] + b1 + bar[i+1:j] + b2 +bar[j+1:]
                            barcodes_mis_dict[bar_mis] = (hamming2(bar_mis, bar), bar)
        else:
            print('mismatch should be less than 3 when matching 8bp barcodes!')
            exit()
    return barcodes_mis_dict

def hamming2(s1, s2):
    """Calculate the Hamming distance between two bit strings
    from https://stackoverflow.com/questions/31007054/hamming-distance-between-two-binary-strings-not-working
    """
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

@click.command(options_metavar='-i METADATA -o DIR [-b barcode] [-m mismatch] [-p THREAD] [-f]',
    short_help='extractor')
@click.option('-i','--metadata', metavar='METADATA', nargs=1, required=True,
    help = 'a tab-delimited txt file including 4 columns of sampleName, fastq1, fastq2 and barcodes (comma-delimited) to extract')
@click.option('-o','--outputdir', default='./extractor', show_default=True,
    metavar='OUTPUTDIR', nargs=1, required=False, 
    type=click.Path(exists=False, file_okay=False, resolve_path=True, writable=True, readable=True),
    help = 'output folder; existed folders is not allowed.')
@click.option('-b','--barcode', #default='96-8bp-barcode', show_default=True,
    metavar='BARCODE', nargs=1, required=False,
    help = 'Use the built-in file if not specify. A file including all available barcodes (each per line). ')
@click.option('-m','--mismatch', default=1, show_default=True, type = int, 
    metavar='MISMATCH', nargs=1, required=False,
    help = 'max number of bases when barcode match')
@click.option('-p','--thread', default=1, show_default=True, type = int, 
    metavar='THREAD', nargs=1, required=False,
    help = 'run in parallel-mode if THREAD > 1')
@click.option('-f','--force', required=False, is_flag=True,
    help = 'Do NOT use it unless you know what you are doing. You may firstly delete the OUTPUT dir if exists. ')
@click.version_option()


def main(metadata, outputdir, barcode, mismatch, thread, force):
    """
    Extract fastq by barcodes for tag-based scRNA-Seq (from TangLab) data.
    version 0.1 (multiprocessing is supported.)
    
    Dependency(including but not limited to):
       Biopython
    """
    logging.basicConfig(filename = './extractor'+time.strftime("%Y-%m-%d-%H_%M_%S",time.localtime(time.time())) +'.log', 
        level = logging.INFO, filemode = 'a', format = '%(asctime)s - %(levelname)s: %(message)s')

    # parse metadata
    fastqList = {}
    with open(metadata, 'r') as meta:
        for i in meta:
            fastqList[i] = i.strip().split('\t')
    
    # get barcodes
    if not barcode:
        run_path = os.path.split(os.path.realpath(__file__))[0]
        barcode = os.path.join(run_path, "96-8bp-barcode")
    barcodes = [s.strip() for s in open(barcode, 'r').readlines()]
    
    # output dir
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    else:
        if(force):
            click.echo('Override the exsited dir!\n')
        else:
            click.echo('output dir already exsited!\n')
            exit()

    # check and assign threads 
    thread = len(fastqList) if thread > len(fastqList) else thread
    if (thread > 0.8* multiprocessing.cpu_count()):
        click.echo('THREAD exceed 80% maximum core numbers! Number '+ str(int(0.5*int( multiprocessing.cpu_count()))) + ' is recommended if available!')
        exit()
    else:
        click.echo('{} threads will be used.'.format(thread))
    ############################
    logging.info('Samples to be processed:\n' + ''.join(fastqList.keys()))

    #using multiprocessing to remedy time-consuming fastq barcode-spliting and htseq-count, which can only be run in single thread mode and typically takes 5 hours one sample.
    mp = {}
    #pairedFastqExtractor(fq1, fq2, fq, outdir, barcodes, mismatch, logging)
    for s, p in fastqList.items():
        while(len(multiprocessing.active_children()) > thread-1):
            time.sleep(10)
        mp[s] = multiprocessing.Process(target=pairedFastqExtractor, args=(p[1].split(','), p[2].split(','), p[0], outputdir, [barcodes[int(i)-1] for i in p[3].split(',')], mismatch, logging))
        mp[s].start()
    for m in mp.values():
        m.join()
    click.echo('All DONE.')
    logging.info('All DONE.')


if __name__ == '__main__':
    main()