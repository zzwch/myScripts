# CentOS镜像使用帮助
http://mirrors.163.com/.help/centos.html

# Download Shiny Server for Red Hat/CentOS 6+
https://rstudio.com/products/shiny/download-server/redhat-centos/

# Shiny Server Professional v1.5.15 Administrator's Guide
https://docs.rstudio.com/shiny-server/   
sudo systemctl start shiny-server   
sudo systemctl status shiny-server   

# Configuration
vi /etc/shiny-server/shiny-server.conf   

# Firewal 防火墙设置 (Centos 7) 
https://blog.csdn.net/f820306455/article/details/84390444   
firewall-cmd --state   
firewall-cmd --list-port   
firewall-cmd --add-port=8080/tcp --permanent #永久开放端口   
firewall-cmd --reload    

