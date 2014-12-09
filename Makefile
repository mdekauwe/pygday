# Mainly because I can't be bothered typing python setup blah
build:
	python setup.py version
	python setup.py build
	python setup.py test
install:
	cat files.txt | xargs rm -rf
	python setup.py version
	python setup.py install --record files.txt
	python setup.py test
clean:
	python setup.py clean 
	cat files.txt | xargs rm -rf
	
