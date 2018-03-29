import click
import os, re  
from multiprocessing import cpu_count
import logging

@click.command(options_metavar='-i INPUT -o OUTPUT -g INDEX [-s ANALYSIS_STEP] [-m] [-b] [-p]',
    short_help='chipliu')
@click.option('-i','--input', metavar='INPUT', nargs=1, required=True, 
    type=click.Path(exists=True, file_okay=False, resolve_path=True),
    help='input data folder; Must contains a folder per sample with input files')
@click.option('-o','--output', metavar='OUTPUT', nargs=1, required=True,
    type=click.Path(exists=False, file_okay=False, resolve_path=True),
    help = 'output folder')
@click.option('-g','--genome', metavar='INDEX', nargs=1, required=True,
    help = 'bowtie2 index file name used for mapping')
@click.option('-s','--steps', metavar='ANALYSIS_STEP', 
    help = 'run only a subset of the workflow; if not specified the complete workflow is run')
@click.option('-m','--macs2_args', metavar='MACS2_ARGS', nargs=1, required=True, default = '-g mm', show_default=True,
    help = 'customize arguments/options of macs2 (peak_calling step). MUST be quoted[IMPORTANT]! Eg. "-g mm --keep_dup auto --pvalue 1e-5 --borad".')
@click.option('-b','--bamCoverage_args', metavar='bamCoverage_ARGS', nargs=1, required=True, default='--normalizeUsing RPGC --effectiveGenomeSize 2494787188 --ignoreForNormalization chrX --extendReads', show_default=True,
    help = 'customize arguments/options of bamCoverage (normalization step). MUST be quoted[IMPORTANT]!')
@click.option('-p','--parallel', metavar='N_CPUS', nargs=1, 
    type=click.INT, default = 1, show_default = True,
    help = 'if specified run ChIP on the parallel-mode')
@click.version_option()
def steps(input, output, genome, steps, macs2_args, bamcoverage_args, parallel):
    """
    A simple command line tool for chipseq data analysis.
    Support for paired-end, illumina 1.9+ phred33
    
	\b
    NOTICE: The current version does not support selective steps.
            The following 5 steps will be executed sequentially.
    ANALYSIS_STEPS:
                 qc: trim apapters and reads quality control
                      - trim_galore + fastqc, 
                      - [REQUIRED] fastq files
            mapping: perform reads alignment
                      - bowtie2 + samtools,
                      - [REQUIRED] clean fastq files
       post_mapping: manipulate mapped bam files
                      - samtools sort + view + rmdup
                      - get uniq mapped reads, and remove duplicates
                      - [REQUIRED] mapped bam files
       peak_calling: identify potential binding sites
                      - macs2
                      - [REQUIRED] deduplicated unique mapped bam files
      normalization: get bigwig files for visualization
                      - deeptools bamCoverage
                      - [REQUIRED] deduplicated unique mapped bam files
    """
    if os.path.isdir(output):
        click.echo('OUTPUT directory exists already! May you wanna try to resume it, otherwise remove/delete it mannually firstly!')
        exit()

    os.makedirs(output)
    flog = open(os.path.join(output, 'chipliu.log'), "w")
    flog.write("A simple command line tool for chipseq data analysis.\n")
    flog.close()
 
    logging.basicConfig(filename = os.path.join(output, 'chipliu.log'), 
        level = logging.INFO, filemode = 'a', format = '%(asctime)s - %(levelname)s: %(message)s')
    #logging.debug('debug')
    #logging.info('info')
    #logging.warning('warn')
    #logging.error('error')

    logging.info('STEP0 - check parameters')
    # check_options(input, output, genome, steps, parallel)
    # check input paired-end fastaq
    logging.info(' '*2 + 'check INPUT: paired-end fastq files for each sample')
    logging.info(' '*4 + input)
    samples = os.listdir(input)
    #click.echo(samples)
    sampleinfos = {}
    re_fastq = re.compile('^(.*)[-_\.]r?(read)?[12]\.f.*q(\.gz)?$',re.I)
    re_fastq1 = re.compile('^(.*)[-_\.]r?(read)?1\.f.*q(\.gz)?$',re.I)
    re_fastq2 = re.compile('^(.*)[-_\.]r?(read)?2\.f.*q(\.gz)?$',re.I)
    for s in samples:
        s_dir = os.path.join(input,s)
        if not os.path.isdir(s_dir):
            continue
        
        #logging.info(s_dir)
        s_files = os.listdir(s_dir)
        s_fastq1 = {}
        s_fastq2 = {}
        for s_file in s_files:
            if not os.path.isfile(os.path.join(s_dir,s_file)):
                continue
            s_re1 = re_fastq1.match(s_file)
            s_re2 = re_fastq2.match(s_file)
            if s_re1 != None:
                s_fastq1[s_re1.group(1)] = s_file
            elif s_re2 != None:
                s_fastq2[s_re2.group(1)] = s_file
        logging.info(' '*6 + s)
        for s_fastq in list(set(s_fastq1.keys()).intersection(set(s_fastq2.keys()))):
            click.echo(sampleinfos)
            if (s in sampleinfos.keys()):
                sampleinfos[s][0].append(s_fastq1[s_fastq])
                sampleinfos[s][1].append(s_fastq2[s_fastq])
            else:
                sampleinfos[s] = [[s_fastq1[s_fastq]], [s_fastq2[s_fastq]]]
            logging.info(' '*8 + s_fastq1[s_fastq])
            logging.info(' '*8 + s_fastq2[s_fastq])
        click.echo(sampleinfos)
    # check output existence
    logging.info(' '*2 + 'check OUTPUT: existence')
    if os.path.isdir(output):
        logging.info(' '*4 + output + ': OUTPUT directory has been created!')
    else:
        logging.info(' '*4 + ': OUTPUT directory cannot be created!?')
        exit()
    # check genome existence
    logging.info(' '*2 + 'check INDEX: existence')
    if (os.path.isfile(genome+'.1.bt2') and os.path.isfile(genome+'.2.bt2') and os.path.isfile(genome+'.3.bt2') and os.path.isfile(genome+'.4.bt2')):
        logging.info(' '*4 + 'genome bowtie2 INDEX files exists!')
    else:
        logging.info(' '*4 + 'genome bowtie2 INDEX files do not exist! Please build them first!')
        #exit()
    # check steps 
    logging.info(' '*2 + 'check ANALYSIS_STEP: legality')
    logging.info(' '*4 + 'these analysis will be performed sequentially: qc, mapping, post_mapping, peak_calling, normalization')
    # check parallel
    logging.info(' '*2 + 'check N_CPUS: maximum')
    if (parallel > cpu_count()):
        logging.info(' '*4 + 'N_CPUS exceed maximum core numbers! Number '+ str(int(0.8*int(cpu_count()))) + ' is recommended if available!')
        exit()
    else:
        logging.info(' '*4 + str(parallel) + ' cores will be used for next analyses.')
    
    
    logging.info('STEP0 - END')
    
    ############################
    psamples = sampleinfos.keys()
    logging.info('Samples to be processed:\n' + ', '.join(psamples))
    os.symlink(input, os.path.join(output, 'rawdata'))
    
    clean_dir = os.path.join(output, 'clean_data')
    qc_dir = os.path.join(output, 'clean_fastqc')
    mapping_dir = os.path.join(output, 'mapping_results')
    peak_dir = os.path.join(output, 'peak_results')
    bigwig_dir = os.path.join(output, 'bigwig_results')
    
    os.mkdir(clean_dir)
    os.mkdir(qc_dir)
    os.mkdir(mapping_dir)
    os.mkdir(peak_dir)
    os.mkdir(bigwig_dir)
    
    for p in psamples:
        p_input_dir = os.path.join(input,p)
        p_clean_dir = os.path.join(clean_dir, p)
        p_qc_dir = os.path.join(qc_dir, p)
        p_mapping_dir = os.path.join(mapping_dir, p)
        p_peak_dir = os.path.join(peak_dir, p)
        p_bigwig_dir = os.path.join(bigwig_dir, p)
        
        logging.info(p + ' >>> proceesing start...')
        
        #trim_galore
        os.mkdir(p_clean_dir)
        os.mkdir(p_qc_dir)
        logging.info(p+' >>> qc: trim_galore - cleaning data (trim adapter and remove low quality) and fastqc')
        # click.echo(sampleinfos)
        # click.echo(sampleinfos[p])
        p_fq = [os.path.join(p_input_dir,rv) for r in zip(sampleinfos[p][0],sampleinfos[p][1]) for rv in r]
        p_cmd_qc = 'trim_galore --phred33 --fastqc --fastqc_args "-t %d -o %s" --stringency 10 --length 30 --max_n 15 --trim-n -o %s --paired %s' %(parallel, p_qc_dir, p_clean_dir, ' '.join(p_fq))
        p_ret = os.system(p_cmd_qc)
        if (p_ret != 0):
            logging.warning(p+' is skipped, '+'due to nonzero return!')
            continue
        
        #bowtie2
        os.mkdir(p_mapping_dir)
        logging.info(p+' >>> mapping: bowtie2 - clean reads mapping to genome')
        re_fq_val_1 = re.compile('^(.*)[-_\.]r?(read)?1_val_1\.f.*q(\.gz)?$',re.I)
        re_fq_val_2 = re.compile('^(.*)[-_\.]r?(read)?2_val_2\.f.*q(\.gz)?$',re.I)
        p_fq_val_1 = [os.path.join(p_clean_dir, r) for r in os.listdir(p_clean_dir) if re_fq_val_1.match(r)]
        p_fq_val_2 = [os.path.join(p_clean_dir, r) for r in os.listdir(p_clean_dir) if re_fq_val_2.match(r)]
        p_cmd_bowtie2 = 'bowtie2 -p %s -x %s -1 %s -2 %s -S %s > %s 2>&1' %(parallel, genome, ','.join(p_fq_val_1), ','.join(p_fq_val_2), os.path.join(p_mapping_dir, p+'.mapped.sam'), os.path.join(p_mapping_dir, p+'.mapping.log'))
        p_ret = os.system(p_cmd_bowtie2)
        if (p_ret != 0):
            logging.warning(p+' is skipped, '+'due to nonzero return!')
            continue
        click.echo('bowtie2')
        # samtools
        logging.info(p+' >>> mapping: samtools view - sam to bam')
        p_cmd_samtools = 'samtools view -h -b -@ %d -o %s %s' %(parallel-1, os.path.join(p_mapping_dir, p+'.mapped.bam'), os.path.join(p_mapping_dir, p+'.mapped.sam'))
        p_ret = os.system(p_cmd_samtools)
        if (p_ret != 0):
            logging.warning(p+' is skipped, '+'due to nonzero return!')
            continue
        os.remove(os.path.join(p_mapping_dir, p+'.mapped.sam'))
        
        logging.info(p+' >>> post_mapping: samtools sort - sort bam')
        p_cmd_samtools = 'samtools sort -@ %d -o %s %s' %(parallel-1, os.path.join(p_mapping_dir, p+'.mapped.sorted.bam'), os.path.join(p_mapping_dir, p+'.mapped.bam'))
        p_ret = os.system(p_cmd_samtools)
        if (p_ret != 0):
            logging.warning(p+' is skipped, '+'due to nonzero return!')
            continue
        logging.info(p+' >>> post_mapping: samtools view - get unique mapped bam')
        p_cmd_samtools = 'samtools view -h -b -q 30 -@ %d -o %s %s' %(parallel-1, os.path.join(p_mapping_dir, p+'.mapped.sorted.uniq.bam'), os.path.join(p_mapping_dir, p+'.mapped.sorted.bam'))
        p_ret = os.system(p_cmd_samtools)
        if (p_ret != 0):
            logging.warning(p+' is skipped, '+'due to nonzero return!')
            continue
        logging.info(p+' >>> post_mapping: samtools rmdup - remove duplicates')
        p_cmd_samtools = 'samtools rmdup %s %s >> %s' %(os.path.join(p_mapping_dir, p+'.mapped.sorted.uniq.bam'), os.path.join(p_mapping_dir, p+'.mapped.sorted.uniq.dedup.bam'), os.path.join(p_mapping_dir, p+'.mapping.log'))
        p_ret = os.system(p_cmd_samtools)
        if (p_ret != 0):
            logging.warning(p+' is skipped, '+'due to nonzero return!')
            continue
        os.remove(os.path.join(p_mapping_dir, p+'.mapped.bam'))
        os.remove(os.path.join(p_mapping_dir, p+'.mapped.sorted.bam'))
        os.remove(os.path.join(p_mapping_dir, p+'.mapped.sorted.uniq.bam'))
        
        logging.info(p+' >>> post_mapping: samtools index ')
        p_cmd_samtools = 'samtools index -@ %d %s' %(parallel-1, os.path.join(p_mapping_dir, p+'.mapped.sorted.uniq.dedup.bam'))
        p_ret = os.system(p_cmd_samtools)
        if (p_ret != 0):
            logging.warning(p+' is skipped, '+'due to nonzero return!')
            continue
        
        # macs2
        os.mkdir(p_peak_dir)
        logging.info(p+' >>> peak_calling: macs2 callpeak')
        
        p_cmd_macs2 = 'macs2 callpeak -t %s --outdir %s -n %s %s > %s 2>&1' %(os.path.join(p_mapping_dir, p+'.mapped.sorted.uniq.dedup.bam'), p_peak_dir, p, macs2_args, os.path.join(p_peak_dir, p+'.callpeak.log'))
        p_ret = os.system(p_cmd_macs2)
        if (p_ret != 0):
            logging.warning(p+' is skipped, '+'due to nonzero return!')
            continue
        
        # bamCoverage
        os.mkdir(p_bigwig_dir)
        logging.info(p+' >>> normalization: bamCoverage')
        p_cmd_bamCoverage = 'bamCoverage -p %d -b %s -o %s %s > %s 2>&1' %(parallel, os.path.join(p_mapping_dir, p+'.mapped.sorted.uniq.dedup.bam'), os.path.join(p_bigwig_dir, p+'.bigwig'), bamcoverage_args, os.path.join(p_bigwig_dir, p+'.bigwig.log'))
        p_ret = os.system(p_cmd_bamCoverage)
        if (p_ret != 0):
            logging.warning(p+' is skipped, '+'due to nonzero return!')
            continue
        
if __name__ == '__main__':
    steps()
