import click
import os,re
import multiprocessing,time
@click.command(options_metavar='-d DIR -c MD5FILE [-p THREAD] [-r]',
    short_help='md5check')
@click.option('-c','--check', metavar='MD5FILE', nargs=1, required=True, 
    type=click.STRING, default = 'MD5_.*.txt', show_default = True,
    help='A regular expression matching with MD5SUM files.')
@click.option('-d','--dir', metavar='DIR', nargs=1, required=True, 
    type=click.Path(exists = True, resolve_path = True),# show_default = True, default = ".",
    help='check MD5SUM of files under this folder')
@click.option('-p','--thread', metavar='THREAD', nargs=1, required=False,
    type = click.INT, default = 1, show_default = True,
    help = 'run in parallel-mode if THREAD > 1')
@click.option('-r','--recursive', required=False, is_flag=True, show_default = True,
    help = 'search also subfolders')
@click.version_option()
def check(check, dir, thread, recursive):
    """
    Designed for checking md5sum of fastq files from Novogene.
    
    version 0.1 @180914    
    
    \b
    multiprocessing is supported!!! 
    """
    
    sampleinfos = get_samples(dir, check, recursive, {})
    #click.echo(sampleinfos)
    for key,value in sampleinfos.items():
        print('{key}:{value}'.format(key = key, value = value))
    click.echo('Performing md5sum check... this is time consuming')
    mp = {}
    for k,v in sampleinfos.items():
        click.echo(k)
        for f in v:
            click.echo(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())) + '  '+f)
            mp[os.path.join(k,f)] = multiprocessing.Process(target=checker, args=(k, f))
            mp[os.path.join(k,f)].start()
            while(len(multiprocessing.active_children()) > thread-1):
                time.sleep(10)

    for m in mp.values():
        m.join()
    click.echo('All Done!')

def checker(dir, file):
    os.chdir(dir)
    os.system('md5sum -c ' + file)
    
def get_samples(dir, check, recursive, sampleinfos):
    try:
        files = os.listdir(dir)
    except OSError:
        click.echo("Permission denied: "+dir)
        return sampleinfos
    else:
        re_check = re.compile(check,re.I)
    sampleinfos[dir] = []
    for f in files:
        f_path = os.path.join(dir, f)
        if not os.path.isdir(f_path):
            f_re = re_check.match(f)
            if f_re != None:
                sampleinfos[dir].append(f) 
            else:
                continue
        else:
            if recursive:
                sampleinfos.update(get_samples(f_path, check, recursive, {}))
            else:
                continue
    if sampleinfos[dir] == []:
        sampleinfos.pop(dir)
    return sampleinfos
    
if __name__ == '__main__':
    check()
