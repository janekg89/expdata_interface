import os
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
STATICFILES_STORAGE = 'whitenoise.django.GzipManifestStaticFilesStorage'
# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.9/howto/static-files/
STATIC_ROOT = os.path.join(PROJECT_ROOT, 'staticfiles')
STATIC_URL = '/static/'

# Extra places for collectstatic to find static files.
STATICFILES_DIRS = (
    os.path.join(PROJECT_ROOT, 'static'),
)