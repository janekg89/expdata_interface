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
Database can be filled with bibliography data with corresponding figures, tables, datasets
`expdata_interface/db_create.py` : script to fill database

To run the web interface locally use
```
python manage.py runserver 8001
http://127.0.0.1:8001/expdata_interface/
```

Add user for admin interface
```
python manage.py createsuperuser
```



## heroku deployment
The following files are required for heroku
```
Procfile
venv/
```

## todo
- link local bootstrap
- rename model fields
- write some example tests (for db)

## questions
- what is `base_site.html` template and where is the `base.html` ?
