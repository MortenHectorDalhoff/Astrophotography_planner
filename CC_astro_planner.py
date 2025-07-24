# Own import
from lib.util.clear_outside import fetch_weather_image, open_clear_outside
from lib.Observer import Observer

# UI modules
import tkinter as tk
from tkinter import filedialog, ttk
from ttkthemes import ThemedStyle

# Process line
import ctypes
import sys
import os

# Other modules
from PIL import Image, ImageTk
import numpy as np
import webbrowser
import json

def open_plan_file():
    global current_file_path
    global observer
    global location_entry
    global latitude_entry
    global longitude_entry
    global timezone_entry

    print("Opening plan file...")

    file_path = current_file_path = filedialog.askopenfilename(title="Open Observation Plan",
                                                       defaultextension=".json",
                                                       filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
                                                       )
    print(f"Selected file: {file_path}")
    if file_path:
        current_file_path = file_path
        print(f"Current file path set to: {current_file_path}")

        update_hybrid_button()
        print("Hybrid button updated.")

        with open(file_path, "r") as f:
            data = json.load(f)
        obs_data = data.get("Observer", {})
        print(f"Loaded observer data: {obs_data}")

        location_name = obs_data.get("location_name", "")
        latitude = obs_data.get("latitude", None)
        longitude = obs_data.get("longitude", None)
        timezone = obs_data.get("timezone", "")

        # Update observer object
        observer.get_location(location_name, latitude, longitude, timezone)
        print(f"Observer updated: {observer}")

        # Update UI entry fields
        location_entry.config(state="normal")
        location_entry.delete(0, tk.END)
        location_entry.insert(0, location_name)
        print(f"Location entry updated: {location_name}")

        latitude_entry.config(state="normal")
        latitude_entry.delete(0, tk.END)
        latitude_entry.insert(0, f"{latitude:.6f}" if latitude is not None else "")
        latitude_entry.config(state="readonly")
        print(f"Latitude entry updated: {latitude}")

        longitude_entry.config(state="normal")
        longitude_entry.delete(0, tk.END)
        longitude_entry.insert(0, f"{longitude:.6f}" if longitude is not None else "")
        longitude_entry.config(state="readonly")
        print(f"Longitude entry updated: {longitude}")

        timezone_entry.config(state="normal")
        timezone_entry.delete(0, tk.END)
        timezone_entry.insert(0, timezone)
        timezone_entry.config(state="readonly")
        print(f"Timezone entry updated: {timezone}")

        # get weather images
        print("Fetching weather images...")
        forecast_img, logo_img = fetch_weather_image(observer.latitude, observer.longitude)

        forecast_img_tk_new = ImageTk.PhotoImage(forecast_img)
        forecast_img_label.configure(image=forecast_img_tk_new)
        forecast_img_label.image = forecast_img_tk_new
        print("Forecast image updated.")

        logo_img_tk_new = ImageTk.PhotoImage(logo_img)
        logo_img_label.configure(image=logo_img_tk_new)
        logo_img_label.image = logo_img_tk_new
        print("Logo image updated.")

def save_plan_file():
    global current_file_path
    global observer

    print("Saving plan file...")

    if not current_file_path:
        file_path = filedialog.asksaveasfilename(
            title="Save Observation Plan",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        if not file_path:
            print("No file path selected. Aborting save operation.")
            return
        current_file_path = file_path
        print(f"Selected file path: {file_path}")

    data = {
        "Observer": {
            "location_name": getattr(observer, "location_name", ""),
            "latitude": getattr(observer, "latitude", None),
            "longitude": getattr(observer, "longitude", None),
            "timezone": getattr(observer, "timezone", "")
        }
    }

    with open(current_file_path, "w") as f:
        json.dump(data, f, indent=4)
    print(f"Plan file saved to: {current_file_path}")

def new_file():
    global current_file_path
    global observer
    global location_entry
    global latitude_entry
    global longitude_entry
    global timezone_entry
    global forecast_img_label
    global logo_img_label

    print("Resetting interface")
    current_file_path = None
    observer = Observer()

    location_entry.delete(0, tk.END)

    latitude_entry.config(state="normal")
    latitude_entry.delete(0, tk.END)
    latitude_entry.config(state="readonly")

    longitude_entry.config(state="normal")
    longitude_entry.delete(0, tk.END)
    longitude_entry.config(state="readonly")

    timezone_entry.config(state="normal")
    timezone_entry.delete(0, tk.END)
    timezone_entry.config(state="readonly")

    print("updating weather images")
    forecast_img_label.configure(image=blank_forecast_img_tk)
    forecast_img_label.image = blank_forecast_img_tk

    logo_img_label.configure(image=blank_logo_img_tk)
    logo_img_label.image = blank_logo_img_tk

    update_hybrid_button()

def save_as_file():
    global current_file_path

    print("Saving as new file...")

    file_path = filedialog.asksaveasfilename(
        title="Save Observation Plan",
        defaultextension=".json",
        filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
    )
    if not file_path:
        print("No file path selected. Aborting save operation.")
        return
    current_file_path = file_path
    print(f"Selected file path: {file_path}")
    save_plan_file()
    update_hybrid_button()

def update_hybrid_button():
    global current_file_path

    if current_file_path:
        filename = current_file_path.split("/")[-1]
        hybrid_button.config(text=f"Save [{filename}]", command=save_plan_file)
    else:
        hybrid_button.config(text="Open Plan File", command=open_plan_file)

def set_hand_cursor(event):
    event.widget.config(cursor="hand2")

def set_default_cursor(event):
    event.widget.config(cursor="")

def on_search():
    location = location_entry.get().strip()

    if not location:
        print("No location entered.")
        return

    print(f"Searching for location: {location})")

    try:
        observer.get_location(location_name=location)
        print(f"Observer resolved: {observer}")

        # Update entries with resolved values
        latitude_entry.config(state="normal")
        latitude_entry.delete(0, tk.END)
        latitude_entry.insert(0, f"{observer.latitude:.6f}")
        latitude_entry.config(state="readonly")
        print(f"Latitude entry updated: {observer.latitude}")

        longitude_entry.config(state="normal")
        longitude_entry.delete(0, tk.END)
        longitude_entry.insert(0, f"{observer.longitude:.6f}")
        longitude_entry.config(state="readonly")
        print(f"Longitude entry updated: {observer.longitude}")

        timezone_entry.config(state="normal")
        timezone_entry.delete(0, tk.END)
        timezone_entry.insert(0, observer.timezone)
        timezone_entry.config(state="readonly")
        print(f"Timezone entry updated: {observer.timezone}")


        print("Getting weather images...")
        # get weather images
        forecast_img, logo_img = fetch_weather_image(observer.latitude, observer.longitude)
        print("Weather images fetched.")

        forecast_img_tk_new = ImageTk.PhotoImage(forecast_img)
        forecast_img_label.configure(image=forecast_img_tk_new)
        forecast_img_label.image = forecast_img_tk_new
        print("Forecast image updated.")

        logo_img_tk_new = ImageTk.PhotoImage(logo_img)
        logo_img_label.configure(image=logo_img_tk_new)
        logo_img_label.image = logo_img_tk_new
        print("Logo image updated.")

    except Exception as e:
        print(f"Error: {e}")

def add_target(event):
    # raise the new target frame
    print("Adding new target")
    new_target_frame.config(relief="raised")

    target_info = {
        "name": "Andromeda Galaxy [M31]",
        "ra": "+00h42m44.3s",
        "dec": "+41°16′9″"
    }


    # Create a new frame for the target
    frame = ttk.Frame(targets_frame, style='target.TFrame')

    # Make column 2 expandable
    frame.columnconfigure(2, weight=1)

    def remove_target(event):
        frame.destroy()
        if frame in target_frames:
            target_frames.remove(frame)
        print("Target removed.")

    x_button = ttk.Label(frame, text="✕", style='target_normal.TLabel')
    x_button.grid(row=0, column=2, sticky="ne", padx=2, pady=2)
    x_button.bind("<ButtonRelease-1>", remove_target)

    target_name_label = ttk.Label(frame, text=target_info["name"], style='target_bold.TLabel')
    target_name_label.grid(row=0, column=0, columnspan=2, sticky="w", padx=5, pady=5)

    ra_label = ttk.Label(frame, text=f"RA: {target_info['ra']}", style='target_normal.TLabel')
    ra_label.grid(row=1, column=0, sticky="w", padx=5, pady=5)

    ra_label = ttk.Label(frame, text=f"Dec: {target_info['dec']}", style='target_normal.TLabel')
    ra_label.grid(row=1, column=1, sticky="w", padx=5, pady=5)

    # Insert the new frame above new_target_frame
    index = targets_frame.pack_slaves().index(new_target_frame)
    frame.pack(in_=targets_frame, before=new_target_frame, fill=tk.X, padx=5, pady=2)

    # Keep track of frames if needed
    target_frames.insert(index, frame)

def show_burger_menu(event=None):
    burger_menu.tk_popup(
        menu_frame.winfo_rootx() + burger_button.winfo_x(),
        menu_frame.winfo_rooty() + burger_button.winfo_y() + burger_button.winfo_height()
    )

def show_about():
    about_win = tk.Toplevel(root)
    about_win.iconbitmap(r'lib\icon\icon_36x36.ico')
    about_win.title("About")
    about_win.configure(bg="white")
    about_win.resizable(False, False)

    # Load .ico image and display at top center
    icon_img = Image.open(r'lib\icon\icon_128x128.ico')
    icon_img_tk = ImageTk.PhotoImage(icon_img)
    icon_label = tk.Label(about_win, image=icon_img_tk, bg="white")
    icon_label.image = icon_img_tk  # Keep reference
    icon_label.pack(pady=(10, 0))

    message = (
        "Cosmic Curiosity Astronomical Observation Planner\n"
        "Version 0.2.1 Beta\n\n"
        "Developed by Morten Hector Dalhoff\n"
        "Contact: mhd@down-to-earth-media.com\n"
        "Weather data from clearoutside.com\n\n"
        "License:\n"
        "This software is open source and may be used, copied, \n"
        "and modified for personal, academic or professional \n"
        "purposes. Any modified versions must include attribution \n"
        "to the original author, Morten Hector Dalhoff. Commercial \n"
        "redistribution or sale of the software, whether original \n"
        "or modified, is not permitted without prior written permission.\n\n"
        "© 2025 Morten Hector Dalhoff. All rights reserved."
    )

    tk.Label(about_win, text=message, justify="left", bg="white", font=("Arial", 10)).pack(padx=10, pady=10)

    def open_github(event):
        webbrowser.open("https://github.com/MortenHectorDalhoff/Astrophotography_planner")

    def open_clearoutside(event):
        webbrowser.open("https://clearoutside.com")

    github_label = tk.Label(about_win, text="GitHub Repository", fg="blue", cursor="hand2", bg="white",
                            font=("Arial", 10, "underline"))
    github_label.pack(anchor="w", padx=10)
    github_label.bind("<Button-1>", open_github)

    clearoutside_label = tk.Label(about_win, text="Weather data from clearoutside.com", fg="blue", cursor="hand2",
                                  bg="white", font=("Arial", 10, "underline"))
    clearoutside_label.pack(anchor="w", padx=10)
    clearoutside_label.bind("<Button-1>", open_clearoutside)

def exit_app():
    root.quit()

##############
# Initialize #
##############

# Change process line icon
if sys.platform == "win32":
    myappid = 'down-to-earth-media.CC_astro_planner.MainWindow'  # Arbitrary string
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
    # Optionally, set the icon for the executable if you bundle with PyInstaller

observer = Observer()
current_file_path = None
target_frames = []

###############
# MAIN WINDOW #
###############

# Create window
root = tk.Tk()
root.iconbitmap(r'lib\icon\icon_256x256.ico')
root.configure(background = 'grey14')
style =ThemedStyle(root)
style.set_theme('equilux')

title_font = ("Arial", 16, "bold")

style.configure('Background.TFrame', background='grey14')
style.configure('widget.TFrame', background='grey25',relief="flat")
style.configure('widget.TButton', background='white', foreground='white')
style.configure('widget_header.TLabel', background='grey25', foreground='grey75', font=("Arial", 16, "bold"))
style.configure('widget.TLabel', background='grey25', foreground='grey75', font=("Arial", 10, "normal"))

style.configure('target.TFrame', background='grey35', relief="raised", borderwidth=1)
style.configure('target_normal.TLabel', background='grey35', foreground='grey75', font=("Arial", 10, "normal"))
style.configure('target_bold.TLabel', background='grey35', foreground='grey75', font=("Arial", 12, "bold"))


root.title("Cosmic Curiosity Astronomical Observation Planner")
root.geometry("1200x800")

############
# TOP MENU #
############

# Top menu and file open
menu_frame = ttk.Frame(root, style='Background.TFrame')
menu_frame.pack(fill=tk.X, padx=10, pady=5)

# Add burger menu button
burger_button = ttk.Button(menu_frame, text="☰", width=5, command=show_burger_menu)
burger_button.pack(side=tk.LEFT)

# Create the menu
burger_menu = tk.Menu(root, tearoff=0)
burger_menu.add_command(label="New", command=new_file)
burger_menu.add_command(label="Open", command=open_plan_file)
burger_menu.add_command(label="Save", command=save_plan_file)
burger_menu.add_command(label="Save As", command=save_as_file)
burger_menu.add_separator()
burger_menu.add_command(label="About", command=show_about)
burger_menu.add_command(label="Exit", command=exit_app)

hybrid_button = ttk.Button(menu_frame, text="Open Plan File", command=open_plan_file)
hybrid_button.pack(side=tk.LEFT, padx=10)

############
# OBSERVER #
############

# Observer frame
observer_frame = ttk.Frame(root, style='widget.TFrame')
observer_frame.pack(fill=tk.X, padx=10, pady=5)

# Set row weights: only row 3 expands vertically
observer_frame.rowconfigure(0, weight=0)
observer_frame.rowconfigure(1, weight=0)
observer_frame.rowconfigure(2, weight=0)
observer_frame.rowconfigure(3, weight=0)
observer_frame.rowconfigure(4, weight=1)

# Title
(ttk.Label(observer_frame, text="Observer Information", style='widget_header.TLabel')
    .grid(row=0, column=0, columnspan=5, sticky="w", padx=5, pady=(5, 5)))

# Location entry and search
location_label = ttk.Label(observer_frame, text="Location:", style='widget.TLabel')
location_label.grid(row=1, column=0, sticky="w", padx=5, pady=5)

location_entry = tk.Entry(observer_frame)
location_entry.grid(row=1, column=1,columnspan=3, sticky="ew", padx=2)
location_entry.bind("<Return>", lambda event: on_search())

location_search_button = ttk.Button(observer_frame, text="Search", style ='widget.TButton',)
location_search_button.grid(row=1, column=4, sticky="ew", padx=5)
location_search_button.config(command=on_search)

# weather information
blank_forecast_img = Image.fromarray(np.zeros((145, 672, 4), dtype=np.uint8)).convert("RGBA")
blank_forecast_img_tk = ImageTk.PhotoImage(blank_forecast_img)
forecast_img_label = ttk.Label(observer_frame, style='widget.TLabel',image=blank_forecast_img_tk)
forecast_img_label.grid(row=0, column=6, rowspan=5, padx=5, pady=2, sticky="e")
forecast_img_label.bind(
    "<Button-1>",
    lambda event: open_clear_outside(event, latitude=observer.latitude, longitude=observer.longitude)
)
forecast_img_label.bind("<Enter>", set_hand_cursor)
forecast_img_label.bind("<Leave>", set_default_cursor)

blank_logo_img = Image.fromarray(np.zeros((80, 122, 4), dtype=np.uint8)).convert("RGBA")
blank_logo_img_tk = ImageTk.PhotoImage(blank_logo_img)
logo_img_label = ttk.Label(observer_frame, style='widget.TLabel', image=blank_logo_img_tk)
logo_img_label.grid(row=2, column=4, rowspan=3, padx=2, pady=2, sticky="e")
logo_img_label.bind(
    "<Button-1>",
    lambda event: open_clear_outside(event, latitude=observer.latitude, longitude=observer.longitude)
)
logo_img_label.bind("<Enter>", set_hand_cursor)
logo_img_label.bind("<Leave>", set_default_cursor)

# Latitude and Longitude
latitude_label = ttk.Label(observer_frame, text="Latitude:", style='widget.TLabel')
latitude_label.grid(row=2, column=0, sticky="w", padx=5, pady=10)

latitude_entry = ttk.Entry(observer_frame, width=15)
latitude_entry.grid(row=2, column=1, sticky="w", padx=2, pady=10)
latitude_entry.config(state="readonly")

longitude_label = ttk.Label(observer_frame, text="Longitude:", style='widget.TLabel')
longitude_label.grid(row=2, column=2, sticky="w", padx=5, pady=10)

longitude_entry = ttk.Entry(observer_frame, width=15)
longitude_entry.grid(row=2, column=3, sticky="w", padx=2, pady=10)
longitude_entry.config(state="readonly")

# Timezone
timezone_label = ttk.Label(observer_frame, text="Timezone:", style='widget.TLabel')
timezone_label.grid(row=3, column=0, sticky="w", padx=2, pady=10)
timezone_entry = ttk.Entry(observer_frame,)
timezone_entry.grid(row=3, column=1, columnspan=2, sticky="ew", padx=2, pady=10)
timezone_entry.config(state="readonly")

###########
# TARGETS #
###########

targets_frame = ttk.Frame(root, width=300, style='widget.TFrame')
targets_frame.pack(side=tk.LEFT, padx=(10, 5),  pady=(5, 10), fill=tk.Y)
targets_frame.pack_propagate(False)  # Prevent shrinking to fit contents

ttk.Label(targets_frame, text="Targets", style='widget_header.TLabel').pack(anchor="n", pady=5)

new_target_frame = ttk.Frame(targets_frame, style='target.TFrame')
new_target_frame.pack(fill=tk.X, padx=5, pady=5)
new_target_frame.bind("<Button-1>", lambda event: new_target_frame.config(relief="sunken"))
new_target_frame.bind("<ButtonRelease-1>", add_target)  # Reset text on release

new_target_label = ttk.Label(new_target_frame,
                             text="+",
                             background='grey35',
                             foreground='grey75',
                             font=("Arial", 20, "normal"))
new_target_label.pack(pady=2)


############
# CALENDAR #
############

calendar_frame = ttk.Frame(root, style='widget.TFrame')
calendar_frame.pack(side=tk.LEFT, padx=(5, 10), pady=(5, 10), fill=tk.BOTH, expand=True)

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
