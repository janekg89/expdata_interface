# expdata_interface

Django database interface for experimental datasets consisting of
- publication (pdf)
- tables & figures (png)
- corresponding raw data (csv, excel)


## installation
```
mkvirtualenv expdata_interface
(expdata_interface) $ pip install -r requirements.txt
```

## usage
### How to fill the database with data?
`expdata_interface/db_create.py` : reads data and fills database



## heroku deployment
The following files are required for heroku
```
Procfile
venv/
```

## todo
- link local bootstrap
- rename model fields
- delete obsolete functions db_create

## questions
- what is `base_site.html` template and where is the `base.html` ?
