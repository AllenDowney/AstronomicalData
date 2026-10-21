# copy notebooks
cp soln/0*.ipynb .

# add tags to hide the solutions
python remove_soln.py

#rm -f gd1_data.hdf gd1_isochrone.hdf5 MIST_iso_5fd2532653c27.iso.cmd
#pytest --nbmake 0*.ipynb

git add 0*.ipynb
git commit -m "Updating notebooks"
git push
