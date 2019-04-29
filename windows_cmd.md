# 移动硬盘无法访问“文件或目录损坏且无法读取”，怎样才能修复？
```
chkdsk 盘符: /x /v /f
for example:
chkdsk K: /x /v /f
The parameters for this command are:
/f option will attempt to fix any found errors
/v	Displays the name of each file in every directory as the disk is checked.
/x option will force the volume you’re about to check to be dismounted before the utility begins a scan

chkdsk [<Volume>[[<Path>]<FileName>]] [/f] [/v] [/r] [/x] [/i] [/c] [/l[:<Size>]] [/b]  
see: https://technet.microsoft.com/en-us/library/cc730714(v=ws.11).aspx
```
# 修改路由表分配两张网卡的外网优先级
```
ipconfig
```
> 查看外网网卡的默认网关地址，记录下来
```
route delete 0.0.0.0
```
> 删除两张网卡访问任意IP的默认路由表
```
route -p add 0.0.0.0 mask 0.0.0.0 111.200.102.1
```
> 增加外网网卡网关111.200.102.1，作为外网路由，并存入永久路由表
```
route PRINT
```
> 查看路由表
