from setuptools import setup

with open("README.txt", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='Cellcano',
    version='0.0.1',
    author='Wenjing Ma',
    author_email='wenjing.ma@emory.edu',
    url='https://marvinquiet.github.io/Cellcano/',
    description='Cellcano is for supervised cell type identification (celltyping) in single-cell genomics.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT', 
    packages=['Cellcano'],
    entry_points={
        "console_scripts": ["Cellcano=Cellcano.main:cli"]
        }
)
