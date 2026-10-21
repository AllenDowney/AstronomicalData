# Build the Jupyter book version

# if adding chapters, update _toc.yml

# copy notebooks
cp ../soln/0*.ipynb .

# add tags to hide the solutions
python prep_notebooks.py

# build the HTML version
jb build .

# push it to GitHub
ghp-import -n -p -f _build/html
