# coding utf8
import setuptools
from local_reannotation.versions import get_versions

with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

setuptools.setup(
    name="local_reannotation",
    version=get_versions(),
    author="Yuxing Xu",
    author_email="xuyuxing@mail.kib.ac.cn",
    description="Here's an example of a repository",
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url="https://github.com/SouthernCD/local_reannotation",
    include_package_data = True,

    entry_points={
        "console_scripts": ["reanno = local_reannotation.cli:main"]
    },    

    packages=setuptools.find_packages(),

    install_requires=[
        "yxutil",
        "yxseq",
        "numpy",
    ],

    python_requires='>=3.5',
)