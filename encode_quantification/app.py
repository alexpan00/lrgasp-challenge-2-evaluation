import dash
from flask_caching import Cache
import dash_uploader as du


app = dash.Dash(__name__,
                suppress_callback_exceptions=True)
app.title = 'Quantification'

config = {
    "DEBUG": True,          # some Flask specific configs
    "CACHE_TYPE": "simple",  # Flask-Caching related configs
    "CACHE_DEFAULT_TIMEOUT": 300
}
cache = Cache()
du.configure_upload(app, r"tmp")
cache.init_app(app.server, config=config)
