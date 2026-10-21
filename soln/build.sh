# pip install pytest nbmake

# run nbmake
#pytest --nbmake 0*.ipynb

# push the changes
git add 0*.ipynb
git commit -m "Updating solutions"
git push

# copy the data
DATA="gd1_data.hdf gd1_isochrone.hdf5 MIST_iso_5fd2532653c27.iso.cmd"
cp $DATA ../data
cd ../data
git add $DATA
git commit -m "Updating data"
git push
