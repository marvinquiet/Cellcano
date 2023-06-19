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
        'tensorflow',
        'six',
        'anndata',
        'scanpy',
        'numpy',
        'h5py',
        'keras',
        'rpy2'
        ]
)
