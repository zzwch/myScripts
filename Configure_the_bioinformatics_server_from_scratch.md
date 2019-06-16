# 从头开始配置生物信息学服务器
实验室新到一台戴尔服务器，需要从头开始进行配置和部署工作所需的软件环境。
这是实验室第二台服务器，将会用作生物信息学数据分析之用。
主要的使用场景是对不同组学测序数据的预处理和高级分析。会涉及数据的存储和管理，服务器远程访问，支持多用户的R和Python工作环境（使用Rstudio-server和JupyterHub实现）等。

接下来会结合第一台服务器使用中积累的经验和教训（主要是教训，比如①当时安装的OS版本较老，很多依赖不支持，一些软件安装不了，比如docker；②多用户管理经验不足，root与user的包管理层次没有明确；③不同组学数据分析的工作环境不够独立；④目录结构混乱；等等），进行调整调优，希望在新的服务器能够优化部署和管理。

## 下载文件列表
提前在网络环境较好的场景下，下载好必需的软件包。
### Centos-7
https://www.centos.org/download/
CentOS-7-x86_64-DVD-1810.iso
下载网址 http://isoredirect.centos.org/centos/7/isos/x86_64/
说明：
Linux操作系统可选的发行版较多，常见的如Centos, Ubuntu, etc.
这里推荐选择Centos 7.6 (当前最新版)，可根据喜好自行选择其他发行版。
关于不同版本的异同，自行了解，不赘述。

### Anaconda3 
https://www.anaconda.com/distribution/#download-section
Anaconda3-2019.03-Linux-x86_64.sh
下载链接 https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh
说明：
conda 作为目前最优秀的软件包管理工具，能够方便的进行软件包的安装和管理。通过添加bioconda，conda-forge等软件源，能够快速的获取到常用的生物信息学软件，大大减少在软件包安装上耗费的时间。推荐必装。

### Rstudio Server
https://www.rstudio.com/products/rstudio/download-server/
rstudio-server-rhel-1.2.1335-x86_64.rpm
下载链接 https://download2.rstudio.org/server/centos6/x86_64/rstudio-server-rhel-1.2.1335-x86_64.rpm
说明：
RStudio Server enables you to provide a browser based interface to a version of R running on a remote Linux server, bringing the power and productivity of the RStudio IDE to server-based deployments of R.
安装指令:`sudo yum install rstudio-server-rhel-1.2.1335-x86_64.rpm`

## 刻录OS DVD启动光盘
这里默认使用Windows操作系统，推荐使用Imgburn进行刻录。 
官方网站http://www.imgburn.com/index.php?act=download
下载链接http://download.imgburn.com/SetupImgBurn_2.5.8.0.exe
其他刻录软件亦可。

## 配置网络环境
### 制作网线
https://jingyan.baidu.com/article/b2c186c83eb646c46ff6ff62.html
自行搜索，不赘述。
### 连接服务器
网线一头插在路由器或交换机上（注意，应保证服务器与自己的工作电脑在同一个局域网下，便于远程访问。），另一头插在服务器的网口上。
说明：网线插口处有两个指示灯，绿色常亮代表物理连接有效，黄色常亮代表信号传输正常。如果出现闪烁或不亮，可能是水晶头没做好出现接触不良，或者是网线有问题。网线里的细线一定要捋平顺再插入水晶头，尽量多留出线头，避免细线堆叠造成压线不紧而接触不良。

## 配置服务器
### 配置RAID
两个固态硬盘组RAID1，作为系统盘。
剩下的3块12T硬盘，暂不设置。

### 安装CentOS
安装了GNOME Desktop选项
这里注意配置好网络，方便联网。
### Linux服务器网线直连
因为两台服务器各有一个闲置10Gbps网口， 采用交叉线进行直连。
在其中一台服务器A上手动配置IPv4为：192.168.32.1；255.255.255.0；网关留空；其他留空
另一台服务器B上手动配置IPv4为：192.168.32.2；255.255.255.0；192.168.32.1；其他留空
然后两台服务器上分别连接上网卡，在其中一台服务器的终端上尝试`ping 192.168.32.2` (假设在A上操作)，ping通即为成功。

之后可通过连接服务器（GNOME自带），使用ssh协议，互相进行访问实现文件传输。（由于网线是六类，只能达到1000Mbps的速度，已经算很快了，可以考虑买一根万兆的超六类网线实现1GB/s的传输速度）

### 配置CentOS 7.6 国内软件源
https://mirrors.tuna.tsinghua.edu.cn/help/centos/
https://mirrors.tuna.tsinghua.edu.cn/help/epel/
### 挂载数据硬盘
使用GNOME 3自带的【工具】->【硬盘】，进行硬盘的格式化（ext4），并配置三个硬盘挂载点分别为/data1，/data2, /data3；记得勾选超级权限，避免用户误卸载数据盘

### 安装软件包
#### 安装NTFS支持
需要完成EPEL源的配置
yum install ntfs-3g 
#### 安装exFAT支持
我的U盘是exFAT格式（为了便于在MAC OS上用），这个可以跳过。
可自行百度搜索centos 7 exfat
https://blog.csdn.net/shile/article/details/52202030
#### 安装Anaconda3
以root权限安装
su
输入超级用户密码，然后
bash Anaconda3-2019.03-Linux-x86_64.sh
#### 安装R
因为已经安装了EPEL源，所以可以直接yum安装
yum install R
安装的是3.5.2版本
这里使用yum安装R，不采用`conda install r`，是因为rstudio-server是使用yum安装，若使用conda装会出现一些动态链接库的问题（yum管理的库与conda管理的库不同导致），不想折腾了。
##### 待补充R包安装
从旧服务器上获取package list，然后在新服务器上安装。

#### 安装Rstudio Server
yum install rstudio-server-rhel-1.2.1335-x86_64.rpm
配置/etc/rstudio/rserver.conf
主要是端口
`www-port=47283`

#### 安装Jupyterhub and Jupyterlab
conda install jupyterhub

#### 安装R kernel for Jupyterlab 
pypi 源
https://mirrors.tuna.tsinghua.edu.cn/help/pypi/

#### 安装docker
跟着docker-ce的安装指引走
https://docs.docker.com/install/linux/docker-ce/centos/

安装好，配置下`/etc/docker/daemon.json` (没有的话，新建一个,好像内部文字的引号必须双引号，反正用下面这应该能行）
{
   "data-root": "/data1/root/docker-data/"
}
##### docker pull higlass/higlass-docker
https://github.com/higlass/higlass-docker
##### docker pull epgg/eg-react
https://github.com/lidaof/eg-react
https://eg.readthedocs.io/en/latest/installation.html
### 创建用户
自编脚本批量创建用户及初始化密码，file中每行对应一个用户名。
```
cat file | while read user
do
  useradd -d /data2/$user -m -g Liulab $user
  echo ${user}:123456 > passwd.tmp
  chpasswd < passwd.tmp
done
```
### 安装生产环境
从旧服务器上导出常用软件列表，在新服务器上使用conda进行批量安装。
R 源https://mirrors.tuna.tsinghua.edu.cn/help/CRAN/
-- R-Kernels
