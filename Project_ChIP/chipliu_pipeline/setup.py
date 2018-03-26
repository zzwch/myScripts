from setuptools import setup

setup(
    name='chipliu',
    version='0.1',
    py_modules=['chipliu'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        chipliu=chipliu:steps
    ''',
)