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

## 配置服务器
### 配置RAID

### 安装CentOS
...待补充...

### 配置CentOS 7.6 国内软件源
https://mirrors.tuna.tsinghua.edu.cn/help/centos/
https://mirrors.tuna.tsinghua.edu.cn/help/epel/
### 更新软件包

### 安装Anaconda3
pypi 源
https://mirrors.tuna.tsinghua.edu.cn/help/pypi/
### 安装生产环境
从旧服务器上导出常用软件列表，在新服务器上使用conda进行批量安装。
- R 
R 源https://mirrors.tuna.tsinghua.edu.cn/help/CRAN/
- Rstudio-server
- Jupyterlab
-- R-Kernels
