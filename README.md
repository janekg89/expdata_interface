# expdata_interface

Django database interface for experimental datasets consisting of
- publication (pdf)
- tables & figures (png)
- corresponding raw data (csv, excel)


## Installation
```
mkvirtualenv expdata_interface
(expdata_interface) $ pip install -r requirements.txt
```

## Usage
Create empty database
```
python manage.py makemigrations
python manage.py migrate
```
Add user for admin interface
```
python manage.py createsuperuser
```

Database can be filled with bibliography data with corresponding figures, tables, datasets
`expdata_interface/db_create.py` : script to fill database


To run the web interface locally use
```
python manage.py runserver 8001
http://127.0.0.1:8001/expdata_interface/
```

## heroku deployment
The following files are required for heroku
```
Procfile
venv/
```

## todo
- navigation menu html
- views for figures, tables, csv
- update admin list views
- custom css instead of headings
