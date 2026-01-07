.PHONY: all clean chmod

all: chmod
	$(MAKE) -C src

chmod:
	chmod +x ordopan.py \
	         src/*.py \
	         src/utils/*.py

clean:
	$(MAKE) -C src clean
