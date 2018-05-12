import click
import ConfigParser
import logging
import os, re  
from multiprocessing import cpu_count
from Bio import SeqIO
import gzip
import json

def mymkdir(dir):
    if not os.path.isdir(dir):
        os.mkdir(dir)
def hamming2(s1, s2):
    """Calculate the Hamming distance between two bit strings
    from https://stackoverflow.com/questions/31007054/hamming-distance-between-two-binary-strings-not-working
    """
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))
def bestbarcode(tag, barcodes, mismatch):
    dist = {}
    if tag in barcodes:
        return [0, [tag]]
    if mismatch <= 0:
        return None
    # consider allowing mismatch 
    for b in barcodes:
        d = hamming2(tag, b)
        if d in dist.keys():
            dist[d] = dist[d]+[b]
        else:
            dist[d] = [b]
    d_min = int(min(dist.keys()))
    if d_min > mismatch:
        return None
    return [d_min, dist[d_min]]

def paired2single(fq1, fq2, barcodes, mismatch, fq, others, tso, polya, min_len):
    mismatch = int(mismatch) # in case of str
    tso_n = len(tso)
    out1 = gzip.open(fq, 'w')
    out2 = gzip.open(others, 'w')
    bbcount = dict(zip(barcodes, [[[0 for i in range(mismatch+1)] for j in range(2)] for k in range(len(barcodes))]))
    bbcount['ambiguous'] = [0 for i in range(2)]
    bbcount['unmatched'] = [0 for j in range(2)]
    
    if len(fq1) == len(fq2):
        for i in range(0, len(fq1)):
            in1 = gzip.open(fq1[i],'rU')
            in2 = gzip.open(fq2[i],'rU')
            r1 = SeqIO.parse(in1, 'fastq')
            r2 = SeqIO.parse(in2, 'fastq')
            for r in r2:
                r_flag = 0
                tag = str(r.seq[0:8])
                bb = bestbarcode(tag, barcodes, mismatch)
                if bb == None:
                    bbcount['unmatched'][0] += 1
                    r_flag = -1
                elif len(bb[1]) > 1:
                    bbcount['ambiguous'][0] += 1
                    r_flag = 2
                else:
                    bbcount[bb[1][0]][0][bb[0]] += 1
                    r_flag = 1
                rr = r1.next()
                # trim tso and polya , and skip read less than 50nt
                ind_tso = 0
                ind_polya = len(rr)
                ind_flag = False
                if tso in rr.seq: # using `in` first to save time in case that most reads did not include tso or polya
                    ind_tso = rr.seq.rfind(tso) + tso_n + 3 # sometimes GGG is at the end of tso
                    ind_flag = True
                if polya in rr.seq:
                    ind_polya = rr.seq.find(polya)
                    ind_flag = True
                if ind_flag:
                    if((ind_polya - ind_tso) >= min_len):
                        rr = rr[ind_tso:ind_polya]
                    else:
                        continue
                if r_flag == 1:
                    rr.id = bb[1][0] + str(r.seq[8:16]) +'_'+ str(rr.id)
                    SeqIO.write(rr, out1, 'fastq')
                    bbcount[bb[1][0]][1][bb[0]] += 1
                else:
                    rr.id = str(r.seq[0:16]) + '_' + str(rr.id)
                    SeqIO.write(rr, out2, 'fastq')
                    if r_flag == -1:
                        bbcount['unmatched'][1] += 1
                    elif r_flag == 2:
                        bbcount['ambiguous'][1] += 1
                    else:
                        exit('r_flag is illegal!')
            in1.close()
            in2.close()
    else:
        return None
    out1.close()
    out2.close()
    return bbcount
def click_exit(log):
    click.echo('There is something wrong, check the log file: ' + log)
    exit()
def umi_count(sam, txt, barcodes, ambiguous):
    def parse_tags(tags):
        gene = "__"
        for tag in tags:
            f_tag = tag.split(":")
            if f_tag[0] == "XF" and f_tag[1] == "Z":
                gene = f_tag[2]
        return gene
        
    umimat = {}
    with open(sam,'r') as alignments:
        for a in alignments:
            a = a.strip()
            f = a.split()
            bar = f[0][0:8]
            umi = f[0][8:16]
            tags = f[11:]
            gene = parse_tags(tags)
            if gene[0:2] == "__":
                if (ambiguous and gene[2:11] == 'ambiguous'):
                    gene = gene.split('[')[1].split(']')[0].split('+')
                else:
                    continue
            if ambiguous and type(gene) == list:
                for g in gene:
                    if g not in umimat.keys():
                        umimat[g] = {}
                    if bar not in umimat[g].keys():
                        umimat[g][bar] = {}
                    if umi not in umimat[g][bar].keys():
                        umimat[g][bar][umi] = 0
                    umimat[g][bar][umi] += 1
            else:
                g = gene
                if g not in umimat.keys():
                    umimat[g] = {}
                if bar not in umimat[g].keys():
                    umimat[g][bar] = {}
                if umi not in umimat[g][bar].keys():
                    umimat[g][bar][umi] = 0
                umimat[g][bar][umi] += 1

    genes=umimat.keys()
    with open(txt,"w") as out:
        for gene in sorted(genes):
            line = gene
            for bar in barcodes:
                if bar in umimat[gene].keys():
                    uni = str(len(umimat[gene][bar])) + ':' + str(sum(umimat[gene][bar].values()))
                else:
                    uni = '0:0'
                line = line + '\t' + uni
            out.write(line + '\n')

def get_samples(input, logging):
    # return sampleinfos = {
    #                         sampleA : [
    #                                    [/path/to/sampleA1_1.fq.gz, /path/to/sampleA2_1.fq.gz, ...],
    #                                    [/path/to/sampleA1_2.fq.gz, /path/to/sampleA2_2.fq.gz, ...]
    #                                  ],
    #                         sampleB : [
    #                                    [/path/to/sampleB1_1.fq.gz, /path/to/sampleB2_1.fq.gz, ...],
    #                                    [/path/to/sampleB1_2.fq.gz, /path/to/sampleB2_2.fq.gz, ...]
    #                                  ],
    #                         ...
    #                      }
    logging.info(' '*2 + 'check INPUT: paired-end fastq files for each sample')
    logging.info(' '*4 + input)
    inputs = os.listdir(input)
    sampleinfos = {}
    re_fastq = re.compile('^(.*)[-_\.]r?(read)?[12]\.f.*q(\.gz)?$',re.I)
    re_fastq1 = re.compile('^(.*)[-_\.]r?(read)?1\.f.*q(\.gz)?$',re.I)
    re_fastq2 = re.compile('^(.*)[-_\.]r?(read)?2\.f.*q(\.gz)?$',re.I)
    f_fastq1 = {}
    f_fastq2 = {}
    for s in inputs:
        s_path = os.path.join(input, s)
        if not os.path.isdir(s_path):
            s_re1 = re_fastq1.match(s)
            s_re2 = re_fastq2.match(s)
            if s_re1 != None:
                f_fastq1[s_re1.group(1)] = s
            elif s_re2 != None:
                f_fastq2[s_re2.group(1)] = s
            continue
        
        #logging.info(s_path)
        s_files = os.listdir(s_path)
        s_fastq1 = {}
        s_fastq2 = {}
        for s_file in s_files:
            if not os.path.isfile(os.path.join(s_path,s_file)):
                continue
            s_re1 = re_fastq1.match(s_file)
            s_re2 = re_fastq2.match(s_file)
            if s_re1 != None:
                s_fastq1[s_re1.group(1)] = s_file
            elif s_re2 != None:
                s_fastq2[s_re2.group(1)] = s_file
        logging.info(' '*6 + s)
        for s_fastq in list(set(s_fastq1.keys()).intersection(set(s_fastq2.keys()))):
            if (s in sampleinfos.keys()):
                sampleinfos[s][0].append(os.path.join(s_path, s_fastq1[s_fastq]))
                sampleinfos[s][1].append(os.path.join(s_path, s_fastq2[s_fastq]))
            else:
                sampleinfos[s] = [[os.path.join(s_path, s_fastq1[s_fastq])], [os.path.join(s_path, s_fastq2[s_fastq])]]
            logging.info(' '*8 + s_fastq1[s_fastq])
            logging.info(' '*8 + s_fastq2[s_fastq])
    for f_fastq in list(set(f_fastq1.keys()).intersection(set(f_fastq2.keys()))):
        logging.info(' '*6 + f_fastq)
        if (f_fastq in sampleinfos.keys()):
            click.echo('sample names is duplicated!')
            exit()
            #sampleinfos[f_fastq][0].append(f_fastq1[f_fastq])
            #sampleinfos[f_fastq][1].append(f_fastq2[f_fastq])
        else:
            sampleinfos[f_fastq] = [[os.path.join(input, f_fastq1[f_fastq])], [os.path.join(input, f_fastq2[f_fastq])]]
        logging.info(' '*8 + f_fastq1[f_fastq])
        logging.info(' '*8 + f_fastq2[f_fastq])
    return sampleinfos

def step(s, cmd, run, logging, info):
    #print cmd
    if run:
        ret = os.system(cmd)
        if (ret != 0):
            logging.warning(s +' ::: '+ info +' skipped, due to nonzero return!')
            return None
        else:
            logging.info(s + ' ::: '+info+': finished') # allright
    else:
        logging.info(s +' ::: '+info+': skipped as you wish.') 
    return 'OK'

def config_check(cf):
    #Required SECTIONS
    secs = ['execute_steps', 'options', 'reference', 'annotation'] 
    #Check SECTIONS
    if(len(set(secs).difference(set(cf.sections()))) > 0):
        click.echo('Check your config file. It must contain sections: '+ ', '.join(secs) + '.')
        exit()
    
    #Required OPTIONS
    sec_opt = {'execute_steps' : ['quality_contrl', 'paired2single', 'clean_reads', 'read_mapping', 'gene_quantify', 'unmapped_remapping', 'unmapped_requantify', 'summary_results'],
               'options'       : ['input', 'output', 'sample', 'thread', 'mismatch', 'minlength', 'barcode', 'tso', 'polya', 'adaptor', 'max_n', 'reference', 'gtf', 're_reference', 're_gtf'],
               'reference'     : [],
               'annotation'    : [],
               'tools'         : ['fastqc', 'cutadapt', 'hisat2', 'samtools', 'htseq-count', 'bam2fastx', 'bamtools', "rscript"]}
    #Check OPTIONS
    for sec in sec_opt.keys():
        opt = sec_opt[sec]
        if(sec == 'reference'):
            if(cf.getboolean('execute_steps', 'read_mapping')):
                ref = cf.get('options','reference')
                if ref=='':
                    click.echo('There is no setting for options:reference, check the config file')
                    exit()
                else:
                    opt = opt + [ref]
            if(cf.getboolean('execute_steps', 'unmapped_remapping')):
                ref = cf.get('options','re_reference')
                if ref=='':
                    click.echo('There is no setting for options:re_reference, check the config file')
                    exit()
                else:
                    opt = opt + ref.split(',')
        if(sec == 'annotation'):
            if(cf.getboolean('execute_steps', 'gene_quantify')):
                gtf = cf.get('options','gtf')
                if gtf=='':
                    click.echo('There is no setting for options:gtf, check the config file')
                    exit()
                else:
                    opt = opt + [gtf]
            if(cf.getboolean('execute_steps', 'unmapped_requantify')):
                gtf = cf.get('options','re_gtf')
                if gtf=='':
                    click.echo('There is no setting for options:re_gtf, check the config file')
                    exit()
                else:
                    opt = opt + gtf.split(',')
        sec_opt[sec] = opt
        if(len(set(opt).difference(set(cf.options(sec)))) > 0):
            click.echo('Check your config file. SECTION `' + sec + '` must contain options: '+ ', '.join(opt) + '.')
            exit()
    #check depended tools
    #TO DO
    return sec_opt

@click.command(options_metavar='-c CONFIG [-i INPUT] [-o OUTPUT] [-s SAMPLE] [-p THREAD] [-f]',
    short_help='smartliu')
@click.option('-c','--config', metavar='CONFIG', nargs=1, required=True,
    #type= click.File(mode='r', encoding=None, errors='strict', lazy=None, atomic=False),
    help = 'configuration file or short name, such as mm10, mm10.refgene, hg19, hg19.refgene.')
@click.option('-i','--input', metavar='INPUT', nargs=1, required=False, 
    type=click.Path(exists=True, file_okay=False, resolve_path=True),
    help='input data folder. Must contains paired-end fastq files whose name should match the regular expression of `SAMPLE_[rR]?[12].f(ast)?q(.gz)?` (SAMPLE list is included in file specified by -s option)')
@click.option('-o','--output', metavar='OUTPUT', nargs=1, required=False, 
    type=click.Path(exists=False, file_okay=False, resolve_path=True, writable=True, readable=True),
    help = 'output folder; existed folders is not allowed.')
@click.option('-s','--sample', metavar='SAMPLE', nargs=1, required=False,
    type = click.File(mode='r', encoding=None, errors='strict', lazy=None, atomic=False),
    help = 'use Default to process all matched samples in INPUT or specify a file for the option. Notice, one SAMPLE per line in that specified file, and only included samples will be processed.')
@click.option('-p','--thread', metavar='THREAD', nargs=1, required=False,
    type=click.INT, show_default = False,
    help = 'run in parallel-mode if THREAD > 1')
@click.option('-f','--force', required=False, is_flag=True,
    help = 'delete the OUTPUT dir if exists[caution]. Do NOT use it unless you know what you are doing.')
@click.version_option()


def smart(config, input, sample, output, thread, force):
    """
    A simple command line tool for tag-based scRNA-Seq (from TangLab) data analysis.
    Support for paired-end, illumina 1.9+ phred33

    \b
    Default: 
         input  = ./raw_data/      ; --input, -i     use default or specify it
         output = ./out_smartliu/  ; --output, -o    use default or specify it
         sample =                  ; --sample, -s    use default to un all samples
         thread = 12               ; --thread, -p    use default or specify it
    
    Dependency(including but not limited to):
         fastqc, cutadapt (need python3 for parellel), hisat2, htseq-count, 
       See config file for more configuration information.
    """
    py_path = os.path.split(os.path.realpath(__file__))[0]
    myconfigs = os.listdir(os.path.join(py_path, 'configs'))
    if(config+'.config' in myconfigs):
        config = os.path.join(py_path,'configs', config+'.config')
    cf = ConfigParser.ConfigParser()
    cf.read(config)
    # override config options using command line options
    if(input):
        cf.set('options', 'input', input)
    if(sample):
        cf.set('options', 'sample', sample)
    if(output):
        cf.set('options', 'output', output)
    if(thread):
        cf.set('options', 'thread', thread)
    
    # check output existence
    if os.path.isdir(cf.get('options', 'output')):
        if force:
            click.echo('OUTPUT directory exists already! files in it will be overrided!')
        else:
            click.echo('OUTPUT directory exists already! May you wanna try to resume it, otherwise remove/delete it mannually firstly!')
            exit()
    opts = config_check(cf)
    click.echo('config checked')
    # get option settings 
    cf_step_qc = cf.getboolean('execute_steps', 'quality_contrl')
    cf_step_single = cf.getboolean('execute_steps', 'paired2single')
    cf_step_clean = cf.getboolean('execute_steps', 'clean_reads')
    cf_step_mapping = cf.getboolean('execute_steps', 'read_mapping')
    cf_step_quantify = cf.getboolean('execute_steps', 'gene_quantify')
    cf_step_remapping = cf.getboolean('execute_steps', 'unmapped_remapping')
    cf_step_requantify = cf.getboolean('execute_steps', 'unmapped_requantify')
    cf_step_summary = cf.getboolean('execute_steps', 'summary_results')
    
    cf_opt_input = cf.get('options', 'input')
    cf_opt_output = cf.get('options', 'output')
    cf_opt_sample = cf.get('options', 'sample')
    cf_opt_thread = cf.getint('options', 'thread')
    cf_opt_mismatch = cf.getint('options', 'mismatch')
    cf_opt_minlength = cf.getint('options', 'minlength')
    cf_opt_barcode = cf.get('options', 'barcode')
    cf_opt_tso = cf.get('options', 'tso')
    cf_opt_polya = cf.get('options', 'polya')
    cf_opt_adaptor = cf.get('options', 'adaptor').split(',')
    cf_opt_maxn = cf.getfloat('options', 'max_n')
    cf_opt_ref = cf.get('options', 'reference')
    cf_opt_gtf = cf.get('options', 'gtf')
    cf_opt_reref = cf.get('options', 're_reference').split(',')
    cf_opt_regtf = cf.get('options', 're_gtf').split(',')
    
    cf_tool_fastqc = cf.get('tools', 'fastqc')
    cf_tool_cutadapt = cf.get('tools', 'cutadapt')
    cf_tool_hisat2 = cf.get('tools', 'hisat2')
    cf_tool_samtools = cf.get('tools', 'samtools')
    cf_tool_htseq = cf.get('tools', 'htseq-count')
    cf_tool_bam2fastx = cf.get('tools', 'bam2fastx')
    cf_tool_bamtools = cf.get('tools', 'bamtools')
    cf_tool_rscript = cf.get('tools', 'rscript')
    # just for short
    input = os.path.abspath(cf_opt_input)
    output = os.path.abspath(cf_opt_output)
    mismatch = cf_opt_mismatch
    minlength = cf_opt_minlength
    rscript = cf_tool_rscript
    if cf_opt_barcode == '':
        cf_opt_barcode = os.path.join(py_path, '96-8bp-barcode')
    barcodes = [s.strip() for s in open(cf_opt_barcode, 'r').readlines()]
    if cf_opt_sample != '':
        samples = [s.strip() for s in open(cf_opt_sample, 'r').readlines()]
    # create output dir
    if not os.path.isdir(output):
        os.makedirs(output)
    # back-up the config info for retrospect
    cf.write(open(os.path.join(output, 'smartliu.config'), "w"))
    
    
    # pipeline logging file
    logpath = os.path.join(output, 'smartliu.log')
    logfile = open(logpath, "w")
    logfile.write("A simple command line tool for tag-based scRNA-Seq data analysis.\n")
    logfile.close()
 
    logging.basicConfig(filename = logpath, 
        level = logging.INFO, filemode = 'a', format = '%(asctime)s - %(levelname)s: %(message)s')
    #logging.debug('debug')
    #logging.info('info')
    #logging.warning('warn')
    #logging.error('error')

    logging.info('STEP0 - check parameters')
    # check input paired-end fastaq
    sampleinfos = get_samples(input, logging)
    if cf_opt_sample != '':
        pseudo_sample = set(samples).difference(set(sampleinfos.keys()))
        if(len(pseudo_sample) > 0):
            logging.info('Check your SAMPLE file, inexisitent samples: '+ ', '.join(pseudo_sample) + '.')
            click_exit(logpath)
        else:
            sampleinfos = {key:value for key,value in sampleinfos.items() if key in samples}
    # check steps 
    logging.info(' '*2 + 'check ANALYSIS_STEP: legality')
    cf_steps_bool = [cf_step_qc, cf_step_single, cf_step_clean, cf_step_mapping, cf_step_quantify, cf_step_remapping, cf_step_requantify, cf_step_summary]
    cf_steps = [x for x, y in zip(cf.options('execute_steps'), cf_steps_bool) if y]
    logging.info(' '*4 + 'these analysis will be performed sequentially:' + ', '.join(cf_steps))
    
    # check reference existence
    refs = []
    if cf_step_mapping:
        refs = refs + [cf_opt_ref]
    if cf_step_remapping:
        refs = refs + cf_opt_reref
    if refs != []:
        logging.info(' '*2 + 'check INDEX existence: ' + ', '.join(refs))
        for ref in refs:
            ref_path = cf.get('reference', ref)
            if (os.path.isfile(ref_path+'.1.ht2') and os.path.isfile(ref_path+'.2.ht2') and os.path.isfile(ref_path+'.3.ht2') and os.path.isfile(ref_path+'.4.ht2')):
                continue
            else:
                logging.info(' '*4 + ': reference hisat2 INDEX files do not exist! Please build them first: '+ ref + ' = ' + ref_path)
                click_exit(logpath)
        logging.info(' '*4 + 'reference hisat2 INDEX files exists!')
    # check annotation existence
    gtfs = []
    if cf_step_quantify:
        gtfs = gtfs + [cf_opt_gtf]
    if cf_step_requantify:
        gtfs = gtfs + cf_opt_regtf
    if gtfs != []:
        logging.info(' '*2 + 'check GTF existence: ' + ', '.join(gtfs))
        for gtf in gtfs:
            gtf_path = cf.get('annotation', gtf)
            if (os.path.isfile(gtf_path)):
                continue
            else:
                logging.info(' '*4 + ': annotation GTF files do not exist! Please get them first: '+ gtf + ' = ' + gtf_path)
                click_exit(logpath)
        logging.info(' '*4 + 'annotation GTF files exists!')
    
    # check threads
    logging.info(' '*2 + 'check THREAD: recommended '+ str(0.2*cpu_count()) +', maximum ' + str(cpu_count()))
    if (cf_opt_thread > 0.8*cpu_count()):
        logging.info(' '*4 + 'THREAD exceed 80% maximum core numbers! Number '+ str(int(0.6*int(cpu_count()))) + ' is recommended if available!')
        click_exit(logpath)
    else:
        logging.info(' '*4 + str(cf_opt_thread) + ' cores will be used for next analyses.')
    
    logging.info('STEP0 - END')
    
    ############################
    logging.info('Samples to be processed:\n' + ','.join(sampleinfos.keys()))
    #ret = os.symlink(input, os.path.join(output, 'rawdata'))
    #if(ret !=0 ):
    #    logging.warning('cannot create the link of input in the output dir!')
    #    click_exit(logpath)
    
    outs_dir = os.path.join(output, 'outs')
    outs_fastqc_dir = os.path.join(outs_dir, 'fastqc')
    outs_summary_dir = os.path.join(outs_dir, 'summary')
    
    mymkdir(outs_dir)
    mymkdir(outs_fastqc_dir)
    mymkdir(outs_summary_dir)
    
    smart_dir = output #os.path.join(output, 'smart')
    smart_clean_dir = os.path.join(smart_dir, 'clean_data')
    smart_mapping_dir = os.path.join(smart_dir, 'results_hisat2')
    smart_quantify_dir = os.path.join(smart_dir, 'results_count')

    mymkdir(smart_dir)
    mymkdir(smart_clean_dir)
    mymkdir(smart_mapping_dir)
    mymkdir(smart_quantify_dir)
    
    fq_stat = {} # stat barcodes distribution
    for s,p in sampleinfos.items():
        logging.info(s + ' >>> proceesing start...')
        #0. fastqc
        logging.info(s + ' >>> 0. fastqc: check sequencing quality')
        s_ret = step(s, '%s -t %d -o %s %s' %(cf_tool_fastqc, cf_opt_thread, outs_fastqc_dir, ' '.join(p[0] + p[1])), cf_step_qc, logging, '0. fastqc')
        if s_ret == None:
            click_exit(logpath)
        #1. paird2single
        logging.info(s + ' >>> 1. paired2single ')
        s_fq_valid = os.path.join(smart_clean_dir, s+'_valid.fastq.gz')
        s_fq_others = os.path.join(smart_clean_dir, s+'_others.fastq.gz')
        if cf_step_single:
            s_ret = paired2single(p[0], p[1], barcodes, mismatch, s_fq_valid, s_fq_others, cf_opt_tso, cf_opt_polya, minlength)
            if (s_ret == None):
                logging.warning(s +' ::: 1. paired2single skipped, due to nonzero return!')
                click_exit(logpath)
            else:
                fq_stat[s] = s_ret
                s_barcode_stat = open(os.path.join(smart_clean_dir , s+'.barcodes.json'), 'w')
                s_barcode_stat.write(json.dumps(s_ret))  
                s_barcode_stat.close()  
                logging.info(s + ' ::: 1. paired2single: finished') # allright
        else:
            logging.info(s +' ::: 1. paired2single: skipped as you wish.') 
        #2. cutadapt
        s_fq_valid_trim = os.path.join(smart_clean_dir, s+'_valid_trim.fastq.gz')
        s_log_valid_trim = os.path.join(smart_clean_dir, s+'_valid_trim.log')
        logging.info(s+' >>> 2. cutadapt: cleaning data (trim adapter and remove low quality )')
        s_ret = step(s, '%s -j %s' %(cf_tool_cutadapt, cf_opt_thread) + ' -a '.join(['']+cf_opt_adaptor) + ' --trim-n -m %d --max-n=%d -o %s %s > %s 2>&1' %(minlength, cf_opt_maxn, s_fq_valid_trim, s_fq_valid, s_log_valid_trim), cf_step_clean, logging, '2. cutadapt')
        if s_ret == None:
            click_exit(logpath)
        #3. fastqc
        logging.info(s+' >>> 3. clean-fastqc: check clean quality ')
        s_ret = step(s, '%s -t %d -o %s %s' %(cf_tool_fastqc, cf_opt_thread, outs_fastqc_dir, s_fq_valid_trim), cf_step_qc & cf_step_clean, logging, '3. clean-fastqc')
        if s_ret == None:
            click_exit(logpath)
        #4.1 hisat2
        s_summary_mapping = os.path.join(smart_mapping_dir, s+'.align_summary.txt')
        s_sam_mapping = os.path.join(smart_mapping_dir, s+'.mapped.sam')
        logging.info(s+' >>> 4.1 hisat2: clean reads mapping to reference')
        s_ret = step(s, '%s -p %s --summary-file %s -x %s -S %s -U %s ' %(cf_tool_hisat2, cf_opt_thread, s_summary_mapping, cf.get('reference', cf_opt_ref), s_sam_mapping, s_fq_valid_trim), cf_step_mapping, logging, '4.1 hisat2')
        if (s_ret == None):
            click_exit(logpath)
        # 4.2 samtools view
        s_bam_mapping = os.path.join(smart_mapping_dir, s+'.mapped.bam')
        logging.info(s+' >>> 4.2 samtools view: sam to bam')
        s_ret = step(s, '%s view -h -b --threads %d -o %s %s' %(cf_tool_samtools, cf_opt_thread-1, s_bam_mapping, s_sam_mapping), cf_step_mapping, logging, '4.2 samtools view')
        if (s_ret == None):
            click_exit(logpath)
        if cf_step_mapping:
            os.remove(s_sam_mapping)
        # 4.3 samtools sort
        s_sorted_bam_mapping = os.path.join(smart_mapping_dir, s+'.mapped.sorted.bam')
        logging.info(s+' >>> 4.3 samtools sort: sort bam')
        s_ret = step(s, '%s sort --threads %d -o %s %s' %(cf_tool_samtools, cf_opt_thread-1, s_sorted_bam_mapping, s_bam_mapping), cf_step_mapping, logging, '4.3 samtools view')
        if (s_ret == None):
            click_exit(logpath)
        if cf_step_mapping:
            os.remove(s_bam_mapping)
        #4.4 samtools index
        logging.info(s+' >>> 4.4 samtools index ')
        s_ret = step(s, '%s index -@ %d %s' %(cf_tool_samtools, cf_opt_thread-1, s_sorted_bam_mapping), cf_step_mapping, logging, '4.4 samtools index')
        if (s_ret == None):
            click_exit(logpath)
        
        # 5 htseq-count
        logging.info(s+' >>> 5 htseq-count ')
        s_sam_quantify = os.path.join(smart_quantify_dir, s+'.quantify.sam')
        s_txt_quantify = os.path.join(smart_quantify_dir, s+'.quantify.txt')
        s_ret = step(s, '%s view -h %s | htseq-count -s no -f sam -o%s - %s > %s' %(cf_tool_samtools, s_sorted_bam_mapping, s_sam_quantify, cf.get('annotation', cf_opt_gtf), s_txt_quantify), cf_step_quantify, logging, '5 htseq-count')
        if (s_ret == None):
            click_exit(logpath)
        if cf_step_quantify:
            umi_count(s_sam_quantify, s_txt_quantify, barcodes, False)
            os.remove(s_sam_quantify)
        # 6.1 unmapped reads
        s_fq_unmapped = os.path.join(smart_clean_dir, s+'_valid_trim.unmapped.fastq.gz')
        logging.info(s+' >>> 6.1 unmapped reads')
        s_ret = step(s, '%s -o %s %s' %(cf_tool_bam2fastx, s_fq_unmapped, s_sorted_bam_mapping), cf_step_remapping, logging, '6.1 unmapped reads')
        if (s_ret == None):
            click_exit(logpath)
        if cf_step_remapping or cf_step_requantify:
            for i in range(len(cf_opt_reref)):
                s_reref = cf_opt_reref[i]
                s_reref_path = cf.get('reference', s_reref)
                s_log_remapping = os.path.join(smart_mapping_dir, s+'.unmapped.'+s_reref+'.align_summary.txt')
                s_sam_remapping = os.path.join(smart_mapping_dir, s+'.unmapped.'+s_reref+'.sam')
                s_bam_remapping = os.path.join(smart_mapping_dir, s+'.unmapped.'+s_reref+'.bam')
                s_sorted_bam_remapping = os.path.join(smart_mapping_dir, s+'.unmapped.'+s_reref+'.sorted.bam')
                s_sorted_filtered_bam_remapping = os.path.join(smart_mapping_dir, s+'.unmapped.'+s_reref+'.sorted_filtered.bam')
                s_remapping_cmd = ['' for j in range(6)]
                s_remapping_cmd[0] = '%s -p %d --no-softclip --summary-file %s -x %s --no-unal -S %s -U %s' %(cf_tool_hisat2, cf_opt_thread, s_log_remapping, s_reref_path, s_sam_remapping, s_fq_unmapped)
                s_remapping_cmd[1] = '%s view --threads %d -b -o %s %s' %(cf_tool_samtools, cf_opt_thread, s_bam_remapping, s_sam_remapping)
                s_remapping_cmd[2] = '%s sort --threads %d -o %s %s' %(cf_tool_samtools, cf_opt_thread, s_sorted_bam_remapping, s_bam_remapping)
                s_remapping_cmd[3] = '%s filter -tag XM:0 -in %s -out %s' %(cf_tool_bamtools, s_sorted_bam_remapping, s_sorted_filtered_bam_remapping)
                s_remapping_cmd[4] = '%s index -@ %d %s' %(cf_tool_samtools, cf_opt_thread, s_sorted_filtered_bam_remapping)
                s_remapping_cmd[5] = 'rm %s %s %s' %(s_sam_remapping, s_bam_remapping, s_sorted_bam_remapping)
                
                logging.info(s+' >>> 6.2 unmapped remapping to '+s_reref)
                s_ret = step(s, ' && '.join(s_remapping_cmd), cf_step_remapping, logging, '6.2 unmapped remapping to' + s_reref)
                if s_ret == None:
                    click_exit(logpath)
                s_regtf = cf_opt_regtf[i]
                s_sam_requantify = os.path.join(smart_quantify_dir, s+'.unmapped.'+s_regtf +'.requantify.sam')
                s_txt_requantify = os.path.join(smart_quantify_dir, s+'.unmapped.'+s_regtf +'.requantify.txt')
                logging.info(s+' >>> 6.3 requantify by '+s_regtf)
                s_ret = step(s, '%s view -h %s | %s -s no --nonunique all -f sam -o%s - %s > %s' %(cf_tool_samtools, s_sorted_filtered_bam_remapping, cf_tool_htseq, s_sam_requantify, cf.get('annotation', s_regtf), s_txt_requantify), cf_step_requantify, logging, '6.3 requantify by '+ s_regtf)
                if s_ret == None:
                    logging.warning(s +' ::: skip requantify by ' + s_regtf + ', due to None matched reads?')
                    #click_exit(logpath)
                if cf_step_requantify:
                    umi_count(s_sam_requantify, s_txt_requantify, barcodes, True)
                    os.remove(s_sam_requantify)
    # summary results
    if cf_step_summary:
        s = "Summary"
        logging.info(s+' >>> 7 summary results ')
        rscript_barcode = os.path.join(py_path, 'rscripts','stat_barcodes.R')
        rscript_mapping = os.path.join(py_path, 'rscripts','stat_mapping.R')
        rscript_quantify = os.path.join(py_path, 'rscripts','stat_quantify.R')
        rscript_remapping = os.path.join(py_path, 'rscripts','stat_remapping.R')
        rscript_requantify = os.path.join(py_path, 'rscripts','stat_requantify.R')
        logging.info(s+' >>> 7.1 stat barcodes')
        ret = step(s,'%s %s %s %s %s %s' %(rscript, rscript_barcode, outs_summary_dir, smart_clean_dir, ','.join(sampleinfos.keys()), ','.join(barcodes)), True, logging, '7.1 stat barcodes')
        if ret ==None:
            logging.warning(s +' ::: skip stat_barcodes due to non-zero return?')
            #click_exit(logpath)
        logging.info(s+' >>> 7.2 stat mapping')
        ret = step(s,'%s %s %s %s %s' %(rscript, rscript_mapping, outs_summary_dir, smart_mapping_dir, ','.join(sampleinfos.keys())), True, logging, '7.2 stat mapping')
        if ret ==None:
            logging.warning(s +' ::: skip stat_mapping due to non-zero return?')
            #click_exit(logpath)
        logging.info(s+' >>> 7.3 stat quantify')
        ret = step(s,'%s %s %s %s %s %s' %(rscript, rscript_quantify, outs_summary_dir, smart_quantify_dir, ','.join(sampleinfos.keys()), cf.get('annotation', cf_opt_gtf)), True, logging, '7.3 stat quantify')
        if ret ==None:
            logging.warning(s +' ::: skip stat_quantify due to non-zero return?')
            #click_exit(logpath)
        if cf_step_remapping or cf_step_requantify:
            for i in range(len(cf_opt_reref)):
                s_reref = cf_opt_reref[i]
                s_regtf = cf_opt_regtf[i]
                logging.info(s+' >>> 7.4 stat remapping')
                ret = step(s,'%s %s %s %s %s %s' %(rscript, rscript_remapping, outs_summary_dir, smart_mapping_dir, ','.join(sampleinfos.keys()), s_reref), True, logging, '7.4 stat remapping')
                if ret ==None:
                    logging.warning(s +' ::: skip stat_remapping due to non-zero return?')
                    #click_exit(logpath)
                logging.info(s+' >>> 7.5 stat requantify')
                ret = step(s,'%s %s %s %s %s %s %s' %(rscript, rscript_requantify, outs_summary_dir, smart_quantify_dir, ','.join(sampleinfos.keys()), cf.get('annotation', s_regtf), s_regtf), True, logging, '7.5 stat requantify')
                if ret ==None:
                    logging.warning(s +' ::: skip stat_requantify due to non-zero return?')
                    #click_exit(logpath)
    logging.info('All DONE. Cheers!')
if __name__ == '__main__':
    smart()