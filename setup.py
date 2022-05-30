from setuptools import setup, find_packages
from Cellcano import __version__

with open("DESCRIPTION", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='Cellcano',
    version=__version__,
    author='Wenjing Ma',
    author_email='wenjing.ma@emory.edu',
    url='https://github.com/marvinquiet/Cellcano',
    project_urls={
        "Documentation": "https://marvinquiet.github.io/Cellcano/",
        "Bug Tracker": "https://github.com/marvinquiet/Cellcano/issues",
        },
    description='Cellcano is for supervised cell type identification (celltyping) in single-cell genomics.',
    long_description=long_description,
    long_description_content_type="text/plain",
    license='MIT', 
    packages=find_packages(),
    entry_points={
        "console_scripts": ["Cellcano=Cellcano.main:main"]
        },
    python_requires=">=3.8",
    install_requires=[
        'tensorflow == 2.4.1',
        'six~=1.15.0',
        'anndata == 0.7.4',
        'scanpy == 1.8.2',
        'numpy == 1.19.2',
        'h5py == 2.10.0',
        'keras',
        'rpy2'
        ]
)
