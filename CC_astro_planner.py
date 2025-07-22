# Own import
from lib.Observer import Observer

# UI modules
import tkinter as tk
from tkinter import filedialog, ttk
from ttkthemes import ThemedStyle

# For weather image
import requests
from PIL import Image, ImageTk
import io
import webbrowser
import numpy as np

# Sample data
sample_targets = [
    {"name": "M 51", "ra": "1329m", "dec": "−47°12'"},
    {"name": "VEGA", "ra": "18h 7m", "dec": "+33°47'"},
    {"name": "NGC 2024", "ra": "05h 42m", "dec": "−01°52'"}
]

def open_file():
    filedialog.askopenfilename(title="Open Observation Plan")

def set_hand_cursor(event):
    event.widget.config(cursor="hand2")

def set_default_cursor(event):
    event.widget.config(cursor="")

def on_search():
    location = location_entry.get().strip()
    lat_text = latitude_entry.get().strip()
    lon_text = longitude_entry.get().strip()

    print(f"Searching for location: {location}, Latitude: {lat_text}, Longitude: {lon_text}")

    try:
        latitude = float(lat_text) if lat_text else None
        longitude = float(lon_text) if lon_text else None
    except ValueError:
        latitude = None
        longitude = None

    try:
        if latitude is not None and longitude is not None:
            observer.get_location(latitude=latitude, longitude=longitude)
        elif location:
            observer.get_location(location_name=location)
        else:
            raise ValueError("Please provide a location name or both latitude and longitude.")

        # Update entries with resolved values
        latitude_entry.config(state="normal")
        latitude_entry.delete(0, tk.END)
        latitude_entry.insert(0, f"{observer.latitude:.6f}")
        latitude_entry.config(state="readonly")

        longitude_entry.config(state="normal")
        longitude_entry.delete(0, tk.END)
        longitude_entry.insert(0, f"{observer.longitude:.6f}")
        longitude_entry.config(state="readonly")

        timezone_entry.config(state="normal")
        timezone_entry.delete(0, tk.END)
        timezone_entry.insert(0, observer.timezone)
        timezone_entry.config(state="readonly")

        location_entry.delete(0, tk.END)
        if hasattr(observer, "location_name"):
            location_entry.insert(0, observer.location_name)

        # get weather images
        forecast_img, logo_img = fetch_weather_image(observer.latitude, observer.longitude)

        forecast_img_tk_new = ImageTk.PhotoImage(forecast_img)
        forecast_img_label.configure(image=forecast_img_tk_new)
        forecast_img_label.image = forecast_img_tk_new

        logo_img_tk_new = ImageTk.PhotoImage(logo_img)
        logo_img_label.configure(image=logo_img_tk_new)
        logo_img_label.image = logo_img_tk_new

    except Exception as e:
        print(f"Error: {e}")

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

##############
# Initialize #
##############

title_font = ("Arial", 12, "bold")
observer = Observer()

###############
# MAIN WINDOW #
###############

# Create window
root = tk.Tk()
root.configure(background = 'grey14')
style =ThemedStyle(root)
style.set_theme('equilux')

style.configure('Custom.TFrame', background='grey14')
style.configure('Light.TButton', background='white', foreground='white')

root.title("Cosmic Curiosity Astronomical Observation Planner")
root.geometry("1200x800")

############
# TOP MENU #
############

# Top menu and file open
menu_frame = ttk.Frame(root, style='Custom.TFrame')
menu_frame.pack(fill=tk.X, padx=10, pady=5)
(ttk.Button(menu_frame, text="☰").pack(side=tk.LEFT))
ttk.Button(menu_frame, text="Open Plan File", command=open_file).pack(side=tk.LEFT, padx=10)

############
# OBSERVER #
############

# Observer frame
observer_frame = ttk.Frame(root)
observer_frame.pack(fill=tk.X, padx=10, pady=5)

# Set row weights: only row 3 expands vertically
observer_frame.rowconfigure(0, weight=0)
observer_frame.rowconfigure(1, weight=0)
observer_frame.rowconfigure(2, weight=0)
observer_frame.rowconfigure(3, weight=0)
observer_frame.rowconfigure(4, weight=1)

# Title
(ttk.Label(observer_frame, text="Observer Information", font=title_font)
    .grid(row=0, column=0, columnspan=5, sticky="w", pady=(0, 5)))

# Location entry and search
location_label = ttk.Label(observer_frame, text="Location:")
location_label.grid(row=1, column=0, sticky="w", padx=2, pady=2)

location_entry = tk.Entry(observer_frame)
location_entry.grid(row=1, column=1,columnspan=3, sticky="ew", padx=2)
location_entry.bind("<Return>", lambda event: on_search())

location_search_button = ttk.Button(observer_frame, text="Search", style ='Light.TButton',)
location_search_button.grid(row=1, column=4, sticky="ew", padx=5)
location_search_button.config(command=on_search)

# weather information
black_forecast_img = Image.fromarray(np.zeros((145, 672, 4), dtype=np.uint8)).convert("RGBA")
forecast_img_tk = ImageTk.PhotoImage(black_forecast_img)
forecast_img_label = ttk.Label(observer_frame, image=forecast_img_tk)
forecast_img_label.grid(row=0, column=6, rowspan=5, padx=5, pady=2, sticky="e")
forecast_img_label.bind(
    "<Button-1>",
    lambda event: open_clear_outside(event, latitude=observer.latitude, longitude=observer.longitude)
)
forecast_img_label.bind("<Enter>", set_hand_cursor)
forecast_img_label.bind("<Leave>", set_default_cursor)

black_logo_img = Image.fromarray(np.zeros((80, 122, 4), dtype=np.uint8)).convert("RGBA")
logo_img_tk = ImageTk.PhotoImage(black_logo_img)
logo_img_label = ttk.Label(observer_frame, image=logo_img_tk)
logo_img_label.grid(row=2, column=4, rowspan=3, padx=2, pady=2, sticky="e")
logo_img_label.bind(
    "<Button-1>",
    lambda event: open_clear_outside(event, latitude=observer.latitude, longitude=observer.longitude)
)
logo_img_label.bind("<Enter>", set_hand_cursor)
logo_img_label.bind("<Leave>", set_default_cursor)


# Latitude and Longitude
latitude_label = ttk.Label(observer_frame, text="Latitude:")
latitude_label.grid(row=2, column=0, sticky="w", padx=2, pady=10)

latitude_entry = ttk.Entry(observer_frame, width=15)
latitude_entry.grid(row=2, column=1, sticky="w", padx=2, pady=10)
latitude_entry.config(state="readonly")

longitude_label = ttk.Label(observer_frame, text="Longitude:")
longitude_label.grid(row=2, column=2, sticky="w", padx=2, pady=10)

longitude_entry = ttk.Entry(observer_frame, width=15)
longitude_entry.grid(row=2, column=3, sticky="w", padx=2, pady=10)
longitude_entry.config(state="readonly")

# Timezone
timezone_label = ttk.Label(observer_frame, text="Timezone:")
timezone_label.grid(row=3, column=0, sticky="w", padx=2, pady=10)
timezone_entry = ttk.Entry(observer_frame,)
timezone_entry.grid(row=3, column=1, columnspan=2, sticky="ew", padx=2, pady=10)
timezone_entry.config(state="readonly")

# Add target button
##target_frame = ttk.Frame(root)
##target_frame.pack(fill=tk.X, padx=10, pady=10)
##ttk.Button(target_frame, text="Add Target").pack(side=tk.LEFT)

# Display targets
##targets_frame = ttk.Frame(root)
##targets_frame.pack(side=tk.LEFT, padx=10)

##for target in sample_targets:
##    box = tk.LabelFrame(targets_frame, text=target["name"], padx=5, pady=5)
##   box.pack(pady=5, fill=tk.X)
##   ttk.Label(box, text=f"RA: {target['ra']}").pack(anchor="w")
##   ttk.Label(box, text=f"DEC: {target['dec']}").pack(anchor="w")

# Calendar
##calendar_frame = ttk.Frame(root)
##calendar_frame.pack(side=tk.LEFT, padx=10)

##ttk.Label(calendar_frame, text="April 2024", font=("Arial", 14)).pack()

##calendar = ttk.Treeview(calendar_frame, columns=("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"), show='headings', height=6)
##for col in calendar["columns"]:
##    calendar.heading(col, text=col)
##    calendar.column(col, width=40, anchor='center')

##dates = [str(i+1) for i in range(30)]
##rows = [dates[i:i+7] for i in range(0, len(dates), 7)]
##for row in rows:
##    while len(row) < 7:
##        row.append('')
##    calendar.insert('', tk.END, values=row)

##calendar.pack()

root.mainloop()
