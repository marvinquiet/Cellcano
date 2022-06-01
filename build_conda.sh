if [ -d "./conda-bld" ]
then
    echo "Directory ./conda-bld exists. Clear it."
    rm -rf ./conda-bld/*
else
    echo "Directory ./conda-bld NOT exists. 'mkdir' for it!"
    mkdir ./conda-bld
fi

if [ -x "$(command -v mamba)" ]
then
    echo "'mamba' installed. So we will use The Fast Mamba."
    mybuild="mambabuild"
else
    echo "'mamba' NOT installed. We just use 'conda build'"
    mybuild="build"
fi

# execute build recipe
conda ${mybuild} conda-recipe --output-folder ./conda-bld

# convert for all platforms; this can also be done through anaconda upload
# cd conda-bld
# conda convert -p all osx-64/cellcano*.tar.bz2
# cd ..

# after `anaconda login`
# anaconda upload ./conda-bld/osx-64/cellcano*-0.0.2*.tar.bz2 --all
