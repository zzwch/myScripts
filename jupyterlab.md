# 1. 配置matplotlib 画图支持中文显示
## step1: shell 查看CentOS自带的中文字体
```
fc-list :lang=zh
```
> /usr/share/fonts/wqy-zenhei/wqy-zenhei.ttc
## step2: python 查看matplotlib配置文件位置
```
import matplotlib
print(matplotlib.matplotlib_fname())
```
> /home/Public/BioSoft/anaconda2/envs/python36/lib/python3.6/site-packages/matplotlib/mpl-data/matplotlibrc

## step3 复制CentOS字体文件到matplotlib字体目录下
```
cp /usr/share/fonts/wqy-zenhei/wqy-zenhei.ttc /home/Public/BioSoft/anaconda2/envs/python36/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf/
```
## step4 编辑matplotlib目录中的font_manager.py文件，使matplotlib支持ttc字体
```
vi /home/Public/BioSoft/anaconda2/envs/python36/lib/python3.6/site-packages/matplotlib/font_manager.py
```
> get_fontext_synonyms函数中修改return {'ttf': ('ttf', 'otf')为 return {'ttf': ('ttf', 'otf', 'ttc')

# 2. 配置用户python kernel
在用到的python版本下安装kernel   
```
pip install ipykernel
```   
安装kernel 到jupyterlab   
```
python -m ipykernel install --name py_project_name --user

```
