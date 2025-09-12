import requests
from PIL import Image
import io
import webbrowser

def fetch_weather_image(latitude, longitude):
    """
    Fetch the weather image from Clear Outside based on latitude and longitude.
    """
    # cast latitude and longitude to string with 2 decimal places
    latitude = f"{latitude:.2f}"
    longitude = f"{longitude:.2f}"

    img_url = f"https://clearoutside.com/forecast_image_medium/{latitude}/{longitude}/forecast.png"
    response = requests.get(img_url)
    if response.status_code != 200:
        print("Failed to fetch medium weather image.")
        return None

    img = Image.open(io.BytesIO(response.content))
    # Left, Top, Right, Bottom
    forecast_img = img.crop((0, 80, img.width, img.height))
    logo_img = img.crop((550, 0, img.width, 80))

    return forecast_img, logo_img

def open_clear_outside(event, latitude=None, longitude=None):

    # Format lat and lon to 2 decimal places
    latitude = f"{latitude:.2f}"
    longitude = f"{longitude:.2f}"

    # Open the Clear Outside forecast page in a web browser
    url = f"https://clearoutside.com/forecast/{latitude}/{longitude}"
    webbrowser.open(url)