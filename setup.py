from setuptools import setup, find_packages

setup(
    name="Protein Sequencing",
    version="0.0.1",
    description="Comparison of Protein Sequences",
    install_requires=['numpy', 'matplotlib', 'prettytable', 'click'],
    entry_points="""
    [console_scripts]
    proteinseq=prot_seq:main
    """,
    author="Jonathan Taylor",
    author_email="jonathanataylor1998@yahoo.com",
    packages=find_packages(),
)