# usohelper

Astropy-based utilities to prepare observations


## Requirements / installation

```
conda create -n usohelper25 python=3.12
conda activate usohelper25
conda install -c conda-forge matplotlib=3.10 astroplan=0.10
````

## Examples

Display overview of twilight times etc:
```
python usohelper.py -v
```

Create a program with dithered exposures on M 15:
```
python usohelper.py --targetname "M 15" --prog
```

See help for more options:
```
python usohelper.py -h
```



## Links

https://sbpy.org

https://astroplan.readthedocs.io/en/latest/index.html

(Well-known alternative: https://rhodesmill.org/skyfield/api.html)


https://docs.python.org/3/howto/argparse.html#argparse-tutorial



