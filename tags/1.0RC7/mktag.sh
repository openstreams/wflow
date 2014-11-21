TAG="1.0RC6"
echo $TAG
svn cp -m "made tag: $TAG" https://wflow.googlecode.com/svn/trunk/ https://wflow.googlecode.com/svn/tags/$TAG
