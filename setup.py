from setuptools import setup

with open("README.txt", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='Cellcano',
    version='0.0.2',
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
    packages=['Cellcano'],
    entry_points={
        "console_scripts": ["Cellcano=Cellcano.main:cli"]
        },
    python_requires=">=3.8",
    install_requires=[
        'tensorflow == 2.4.1',
        'anndata == 0.7.4',
        'scanpy == 1.8.2',
        'rpy2',
        'keras',
        'matplotlib'
        ]
)
