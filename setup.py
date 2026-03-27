from setuptools import setup, find_packages


setup(
    name="scpositioner",
    version="1.0.0",
    author='Jingyang Qian',
    author_email='qianjingyang@zju.edu.cn',
    url='https://github.com/ZJUFanLab/scPositioner',
    packages=find_packages(),
    install_requires=[],
    package_data={'': ['*.gz']},
    zip_safe=False,
    python_requires='>=3.9',
    license='GNU General Public License v3.0',
)
