# Own import
from lib.util.clear_outside import fetch_weather_image, open_clear_outside
from lib.Observer import Observer
from lib.AstroTarget import AstroTarget

# UI modules
import tkinter as tk
from tkinter import filedialog, ttk
from ttkthemes import ThemedStyle
import calendar
from datetime import datetime, timedelta

# Process line
import ctypes
import sys
import os

# Other modules
from PIL import Image, ImageTk
import numpy as np
import webbrowser
import json


##### WORK IN PROCESS - NOT FINALIZED #####

#### Main interface functions ####

def add_placeholder(entry, placeholder, color="grey"):
    def on_focus_in(event):
        if entry.get() == placeholder:
            entry.delete(0, tk.END)
            entry.config(fg="black")
    def on_focus_out(event):
        if not entry.get():
            entry.insert(0, placeholder)
            entry.config(fg=color)
    entry.insert(0, placeholder)
    entry.config(fg=color)
    entry.bind("<FocusIn>", on_focus_in)
    entry.bind("<FocusOut>", on_focus_out)

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

#### File operations ####

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

#### cursor functions ####

def set_hand_cursor(event):
    event.widget.config(cursor="hand2")

def set_default_cursor(event):
    event.widget.config(cursor="")

#### Observer and Target Search Functions ####

def on_observer_search():
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

        #scaling image
        scale_factor = 1.2  # Change as needed
        new_width = int(forecast_img.width * scale_factor)
        new_height = int(forecast_img.height * scale_factor)
        forecast_img_resized = forecast_img.resize((new_width, new_height), Image.LANCZOS)


        forecast_img_tk_new = ImageTk.PhotoImage(forecast_img_resized)
        forecast_img_label.configure(image=forecast_img_tk_new)
        forecast_img_label.image = forecast_img_tk_new
        print("Forecast image updated.")

        logo_img_tk_new = ImageTk.PhotoImage(logo_img)
        logo_img_label.configure(image=logo_img_tk_new)
        logo_img_label.image = logo_img_tk_new
        print("Logo image updated.")

    except Exception as e:
        print(f"Error: {e}")

def on_target_search():

    global target

    target_name = target_search_entry.get().strip()

    if not target_name:
        print("No target name entered.")
        return

    print(f"Searching for target: {target_name}")

    try:
        # Create a new AstroTarget instance
        target.resolve_target_from_name(target_name)
        print(f"Target resolved: {target}")

        # Update entries with resolved values
        display_name_entry.config(state="normal")
        display_name_entry.delete(0, tk.END)
        display_name_entry.insert(0, target.pretty_name)
        display_name_entry.config(state="readonly")
        print(f"Display name entry updated: {target.pretty_name}")

        ra_entry.config(state="normal")
        ra_entry.delete(0, tk.END)
        ra_entry.insert(0, target.ra_str)
        ra_entry.config(state="readonly")
        print(f"RA entry updated: {target.ra_str}")

        dec_entry.config(state="normal")
        dec_entry.delete(0, tk.END)
        dec_entry.insert(0, target.dec_str)
        dec_entry.config(state="readonly")
        print(f"Dec entry updated: {target.dec_str}")

    except Exception as e:
        print(f"Error: {e}")

#### Calendar Functions ####

def prev_month():
    global current_year, current_month
    if current_month == 1:
        current_month = 12
        current_year -= 1
    else:
        current_month -= 1
    draw_calendar(canvas, current_year, current_month)

def next_month():
    global current_year, current_month
    if current_month == 12:
        current_month = 1
        current_year += 1
    else:
        current_month += 1
    draw_calendar(canvas, current_year, current_month)

def go_to_today():
    now = datetime.now()
    draw_calendar(canvas, now.year, now.month)

def draw_calendar(canvas, year, month):
    # Set the month/year label
    month_name = calendar.month_name[month]
    month_label.config(text=f"{month_name} {year}")

    canvas.delete("all")

    width = canvas.winfo_width()
    height = canvas.winfo_height()
    cal = calendar.Calendar(firstweekday=0)
    month_days = cal.monthdayscalendar(year, month)
    rows = len(month_days)
    header_h = 30 if rows == 0 else max(20, height // (2 * (rows + 1)))
    cell_h = (height - header_h) // rows if rows else 60
    cell_w = width // 7

    # Parse min/max date
    try:
        min_date = datetime.strptime(min_date_entry.get().strip(), "%Y-%m-%d")
    except Exception:
        min_date = None
    try:
        max_date = datetime.strptime(max_date_entry.get().strip(), "%Y-%m-%d")
    except Exception:
        max_date = None

    # Draw weekday headers
    days_abbr = calendar.day_abbr
    for j, day_name in enumerate(days_abbr):
        x0, y0 = j * cell_w, 0
        canvas.create_rectangle(x0, y0, x0 + cell_w, header_h, fill="#444444", outline="#888888")
        canvas.create_text(
            x0 + cell_w // 2, y0 + header_h // 2,
            text=day_name,
            font=("Arial", 12, "bold"),
            fill="#bfbfbf"
        )

    # Draw days (start from row 1)
    for i, week in enumerate(month_days):
        for j, day in enumerate(week):
            x0, y0 = j * cell_w, header_h + i * cell_h
            x1, y1 = x0 + cell_w, y0 + cell_h
            fill_color = "#323232"  # default

            if day != 0:
                this_date = datetime(year, month, day)
                # Highlight if in range
                if min_date and max_date and min_date <= this_date <= max_date:
                    fill_color = "#666666"  # greenish highlight
                elif (day == datetime.now().day and month == datetime.now().month and year == datetime.now().year):
                    fill_color = "#666666"  # today highlight

            canvas.create_rectangle(x0, y0, x1, y1, fill=fill_color, outline="#888888")
            if day != 0:
                canvas.create_text(
                    x0 + 10, y0 + 10,
                    anchor="nw",
                    text=str(day),
                    font=("Arial", 14, "bold"),
                    fill="#bfbfbf"
                )

def on_canvas_resize(event):
    draw_calendar(canvas, now.year, now.month)

def set_start_date(date_str):
    min_date_entry.delete(0, tk.END)
    min_date_entry.insert(0, date_str)

    # Check and update end date if needed
    end_date = max_date_entry.get().strip()
    try:
        new_start = datetime.strptime(date_str, "%Y-%m-%d")
        if end_date:
            end_dt = datetime.strptime(end_date, "%Y-%m-%d")
            if new_start > end_dt:
                max_date_entry.delete(0, tk.END)
                max_date_entry.insert(0, date_str)
    except ValueError:
        pass  # Ignore invalid date format

    draw_calendar(canvas, current_year, current_month)

def set_end_date(date_str):
    max_date_entry.delete(0, tk.END)
    max_date_entry.insert(0, date_str)

    # Check and update start date if needed
    start_date = min_date_entry.get().strip()
    try:
        new_end = datetime.strptime(date_str, "%Y-%m-%d")
        if start_date:
            start_dt = datetime.strptime(start_date, "%Y-%m-%d")
            if new_end < start_dt:
                min_date_entry.delete(0, tk.END)
                min_date_entry.insert(0, date_str)
    except ValueError:
        pass  # Ignore invalid date format

    draw_calendar(canvas, current_year, current_month)

def on_date_entry_update(event=None):
    draw_calendar(canvas, current_year, current_month)

def on_calendar_right_click(event):
    # Find which day was clicked
    x, y = event.x, event.y
    width = canvas.winfo_width()
    height = canvas.winfo_height()
    cal = calendar.Calendar(firstweekday=0)
    month_days = cal.monthdayscalendar(current_year, current_month)
    rows = len(month_days)
    header_h = 30 if rows == 0 else max(20, height // (2 * (rows + 1)))
    cell_h = (height - header_h) // rows if rows else 60
    cell_w = width // 7

    # Check if click is in the day grid
    if y < header_h:
        return
    row = (y - header_h) // cell_h
    col = x // cell_w
    if 0 <= row < rows and 0 <= col < 7:
        day = month_days[row][col]
        if day == 0:
            return
        date_str = f"{current_year:04d}-{current_month:02d}-{day:02d}"
        # Show context menu
        calendar_menu.delete(0, tk.END)
        calendar_menu.add_command(label="Set Start Date", command=lambda: set_start_date(date_str))
        calendar_menu.add_command(label="Set End Date", command=lambda: set_end_date(date_str))
        calendar_menu.tk_popup(canvas.winfo_rootx() + x, canvas.winfo_rooty() + y)

#### Calculation Functions ####

def on_calculate():
    # 1. Validate observer/location
    if not getattr(observer, "location_name", None) or not getattr(observer, "latitude", None) or not getattr(observer, "longitude", None):
        raise ValueError("Please set a valid observer location.")

    # 2. Validate target
    if not getattr(target, "coordinate", None):
        raise ValueError("Please select a valid target.")

    # 3. Clear all constraints on the target
    if hasattr(target, "clear_constraints"):
        target.clear_constraints()

    # 4. Validate and add constraints
    # Minimum observation hours
    obs_hours = observation_hours_entry.get().strip()
    if obs_hours and obs_hours != "Minimum observation Hours":
        try:
            obs_hours_val = float(obs_hours)
            if obs_hours_val <= 0:
                raise ValueError("Minimum observation hours must be a positive number.")
            target.set_minimum_observation_minutes(int(obs_hours_val * 60))
        except ValueError:
            raise ValueError("Minimum observation hours must be a positive number.")

    # Min altitude
    min_alt = min_altitude_entry.get().strip()
    min_alt_val = None
    if min_alt and min_alt != "0 to 90 (degrees)":
        try:
            min_alt_val = float(min_alt)
            if not (0 <= min_alt_val <= 90):
                raise ValueError("Minimum altitude must be between 0 and 90.")
        except ValueError:
            raise ValueError("Minimum altitude must be a number between 0 and 90.")

    # Max altitude
    max_alt = max_altitude_entry.get().strip()
    max_alt_val = None
    if max_alt and max_alt != "0 to 90 (degrees)":
        try:
            max_alt_val = float(max_alt)
            if not (0 <= max_alt_val <= 90):
                raise ValueError("Maximum altitude must be between 0 and 90.")
        except ValueError:
            raise ValueError("Maximum altitude must be a number between 0 and 90.")

    # Add altitude constraint if valid
    if min_alt_val is not None or max_alt_val is not None:
        vals = []
        if min_alt_val is not None and max_alt_val is not None:
            if min_alt_val > max_alt_val:
                raise ValueError("Minimum altitude cannot be greater than maximum altitude.")
            vals = (min_alt_val, max_alt_val)
        elif min_alt_val is not None:
            vals = (min_alt_val,)
        elif max_alt_val is not None:
            vals = (0, max_alt_val)
        if vals:
            target.add_constraint('altitude', vals)

    # Moon separation
    moon_sep = moon_separation_entry.get().strip()
    if moon_sep and moon_sep != "0 to 180 (degrees)":
        try:
            moon_sep_val = float(moon_sep)
            if not (0 <= moon_sep_val <= 180):
                raise ValueError("Moon separation must be between 0 and 180.")
            target.add_constraint('moon_separation', (moon_sep_val,))
        except ValueError:
            raise ValueError("Moon separation must be a number between 0 and 180.")

    # Moon phase
    moon_phase = moon_phase_entry.get().strip()
    if moon_phase and moon_phase != "0 to 100 (percent)":
        try:
            moon_phase_val = float(moon_phase)
            if not (0 <= moon_phase_val <= 100):
                raise ValueError("Moon phase must be between 0 and 100.")
            target.add_constraint('moon_illumination', (moon_phase_val,))
        except ValueError:
            raise ValueError("Moon phase must be a number between 0 and 100.")

    #5. Validate date range
    min_date_str = min_date_entry.get().strip()
    max_date_str = max_date_entry.get().strip()
    if not min_date_str or not max_date_str:
        # If no dates are set, use the current date
        now = datetime.now()
        min_date_str = now.strftime("%Y-%m-%d")
        max_date_str = now.strftime("%Y-%m-%d")

    try:
        min_date = datetime.strptime(min_date_str, "%Y-%m-%d")
        max_date = datetime.strptime(max_date_str, "%Y-%m-%d")
        if min_date > max_date:
            raise ValueError("Minimum date cannot be after maximum date.")

    except ValueError:
        raise ValueError("Invalid date format. Please use YYYY-MM-DD.")

    # loop through the days in the range
    calendar_events.clear()

    for single_date in (min_date + timedelta(days=n) for n in range((max_date - min_date).days + 1)):
        single_date_observation_window = target.get_observation_window(observer, single_date)
        if single_date_observation_window['is_observable']:
            ovservable_days.append(single_date_observation_window)
            # Extract event info
            start = single_date_observation_window['start']
            end = single_date_observation_window['end']
            total = single_date_observation_window['observable_minutes']
            # Convert astropy Time to datetime if needed
            if hasattr(start, 'datetime'):
                start_dt = start.datetime
            else:
                start_dt = start
            if hasattr(end, 'datetime'):
                end_dt = end.datetime
            else:
                end_dt = end
            # Store event
            calendar_events.append({
                'date': single_date.strftime("%Y-%m-%d"),
                'start': start_dt,
                'end': end_dt,
                'total_minutes': int(total.value) if hasattr(total, 'value') else total
            })
            print(f"Target {target.pretty_name} is observable on {single_date.strftime('%Y-%m-%d')}")

    # Optionally, trigger a calendar redraw here
    draw_calendar(canvas,  , current_month)

##############
# Initialize #
##############

# Change process line icon
if sys.platform == "win32":
    myappid = 'down-to-earth-media.CC_astro_planner.MainWindow'  # Arbitrary string
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
    # Optionally, set the icon for the executable if you bundle with PyInstaller

observer = Observer()
target = AstroTarget()
calendar_events = []
current_file_path = None

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

style.configure("calendar.Treeview", rowheight=75)


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
    .grid(row=0, column=0, columnspan=3, sticky="w", padx=5, pady=(5, 5)))

# Location entry and search
location_label = ttk.Label(observer_frame, text="Location", style='widget.TLabel')
location_label.grid(row=1, column=0, sticky="w", padx=5, pady=2)

location_entry = tk.Entry(observer_frame)
location_entry.grid(row=2, column=0,columnspan=2, sticky="ew", padx=2, pady=(2, 15))
location_entry.bind("<Return>", lambda event: on_observer_search())

location_search_button = ttk.Button(observer_frame, text="Search", style ='widget.TButton',)
location_search_button.grid(row=2, column=2, sticky="nw", padx=5)
location_search_button.config(command=on_observer_search)

# weather information
blank_forecast_img = Image.fromarray(np.zeros((145, 672, 4), dtype=np.uint8)).convert("RGBA")
blank_forecast_img_tk = ImageTk.PhotoImage(blank_forecast_img)
forecast_img_label = ttk.Label(observer_frame, style='widget.TLabel',image=blank_forecast_img_tk)
forecast_img_label.grid(row=1, column=6, rowspan=7, padx=5, pady=2, sticky="nw")
forecast_img_label.bind(
    "<Button-1>",
    lambda event: open_clear_outside(event, latitude=observer.latitude, longitude=observer.longitude)
)
forecast_img_label.bind("<Enter>", set_hand_cursor)
forecast_img_label.bind("<Leave>", set_default_cursor)

blank_logo_img = Image.fromarray(np.zeros((80, 122, 4), dtype=np.uint8)).convert("RGBA")
blank_logo_img_tk = ImageTk.PhotoImage(blank_logo_img)
logo_img_label = ttk.Label(observer_frame, style='widget.TLabel', image=blank_logo_img_tk)
logo_img_label.grid(row=3, column=2, rowspan=4, padx=2, pady=2, sticky="nw")
logo_img_label.bind(
    "<Button-1>",
    lambda event: open_clear_outside(event, latitude=observer.latitude, longitude=observer.longitude)
)
logo_img_label.bind("<Enter>", set_hand_cursor)
logo_img_label.bind("<Leave>", set_default_cursor)

# Latitude and Longitude
latitude_label = ttk.Label(observer_frame, text="Latitude", style='widget.TLabel')
latitude_label.grid(row=3, column=0, sticky="w", padx=5, pady=2)

longitude_label = ttk.Label(observer_frame, text="Longitude", style='widget.TLabel')
longitude_label.grid(row=3, column=1, sticky="w", padx=5, pady=2)

latitude_entry = ttk.Entry(observer_frame, width=15)
latitude_entry.grid(row=4, column=0, sticky="w", padx=2, pady=(2, 15))
latitude_entry.config(state="readonly")

longitude_entry = ttk.Entry(observer_frame, width=15)
longitude_entry.grid(row=4, column=1, sticky="w", padx=2, pady=(2, 15))
longitude_entry.config(state="readonly")

# Timezone
timezone_label = ttk.Label(observer_frame, text="Timezone", style='widget.TLabel')
timezone_label.grid(row=5, column=0, sticky="w", padx=2, pady=2)
timezone_entry = ttk.Entry(observer_frame,)
timezone_entry.grid(row=6, column=0, columnspan=2, sticky="ew", padx=2, pady=(2,15))
timezone_entry.config(state="readonly")

###########
# TARGETS #
###########

targets_frame = ttk.Frame(root, style='widget.TFrame')
targets_frame.columnconfigure(0, weight=1)
targets_frame.columnconfigure(1, weight=1)
targets_frame.pack(side=tk.LEFT, padx=(10, 5),  pady=(5, 10), anchor="nw" , fill=tk.Y)

target_header_label = ttk.Label(targets_frame, text="Target Information", style='widget_header.TLabel')
target_header_label.grid(row=0, column=0, columnspan=4, sticky="w", padx=5, pady=5)

target_search_label = ttk.Label(targets_frame, text="Target Name", style='widget.TLabel')
target_search_label.grid(row=1, column=0, sticky="w", padx=5, pady=2)

target_search_entry = tk.Entry(targets_frame)
target_search_entry.grid(row=2, column=0,columnspan=2, sticky="ew", padx=5, pady=2)
target_search_entry.bind("<Return>", lambda event: on_target_search())

target_search_button = ttk.Button(targets_frame, text="Search", style ='widget.TButton')
target_search_button.grid(row=3, column=0, sticky="ew", padx=5, pady=(2,10))
target_search_button.config(command=on_target_search)

display_name_label = ttk.Label(targets_frame, text="Display Name", style='widget.TLabel')
display_name_label.grid(row=4, column=0, sticky="w", padx=5, pady=2)

display_name_entry = ttk.Entry(targets_frame)
display_name_entry.grid(row=5, column=0, columnspan=2, sticky="ew", padx=5, pady=(2,5))
display_name_entry.config(state="readonly")

ra_label = ttk.Label(targets_frame, text="RA", style='widget.TLabel')
ra_label.grid(row=6, column=0, sticky="ew", padx=5, pady=2)

dec_label = ttk.Label(targets_frame, text="Dec", style='widget.TLabel')
dec_label.grid(row=6, column=1, sticky="ew", padx=5, pady=2)

ra_entry = ttk.Entry(targets_frame, width=15)
ra_entry.grid(row=7, column=0, sticky="ew", padx=5, pady=(2,10))
ra_entry.config(state="readonly")

dec_entry = ttk.Entry(targets_frame, width=15)
dec_entry.grid(row=7, column=1, sticky="ew", padx=5, pady=(2,10))
dec_entry.config(state="readonly")

##############
# CONSTRAINT #
##############

constraint_header_label = ttk.Label(targets_frame, text="Constraints", style='widget_header.TLabel')
constraint_header_label.grid(row=8, column=0, columnspan=2, sticky="w", padx=5, pady=5)

observation_hours_label = ttk.Label(targets_frame, text="Minimum observation Hours", style='widget.TLabel')
observation_hours_label.grid(row=9, column=0, columnspan=2, sticky="w", padx=5, pady=2)
observation_hours_entry = tk.Entry(targets_frame)
observation_hours_entry.grid(row=10, column=0, columnspan=2, sticky="ew", padx=5, pady=(2,5))
add_placeholder(observation_hours_entry, "Minimum observation Hours")

min_altitude_label = ttk.Label(targets_frame, text="Min Altitude", style='widget.TLabel')
min_altitude_label.grid(row=11, column=0, sticky="w", padx=5, pady=2)
max_altitude_label = ttk.Label(targets_frame, text="Max Altitude", style='widget.TLabel')
max_altitude_label.grid(row=11, column=1, sticky="w", padx=5, pady=2)

min_altitude_entry = tk.Entry(targets_frame)
min_altitude_entry.grid(row=12, column=0, sticky="ew", padx=5, pady=(2,5))
add_placeholder(min_altitude_entry, "0 to 90 (degrees)")

max_altitude_entry = tk.Entry(targets_frame)
max_altitude_entry.grid(row=12, column=1, sticky="ew", padx=5, pady=(2,5))
add_placeholder(max_altitude_entry, "0 to 90 (degrees)")

moon_separation_label = ttk.Label(targets_frame, text="Moon Separation", style='widget.TLabel')
moon_separation_label.grid(row=13, column=0, sticky="w", padx=5, pady=2)
moon_separation_entry = tk.Entry(targets_frame)
moon_separation_entry.grid(row=14, column=0, sticky="ew", padx=5, pady=(2,5))
add_placeholder(moon_separation_entry, "0 to 180 (degrees)")

moon_phase_label = ttk.Label(targets_frame, text="Moon Phase", style='widget.TLabel')
moon_phase_label.grid(row=15, column=0, sticky="w", padx=5, pady=0)
moon_phase_entry = tk.Entry(targets_frame)
moon_phase_entry.grid(row=16, column=0, sticky="ew", padx=5, pady=(2,5))
add_placeholder(moon_phase_entry, "0 to 100 (percent)")

############
# CALENDAR #
############

calendar_frame = ttk.Frame(root, style='widget.TFrame')
calendar_frame.pack(side=tk.LEFT, padx=(5, 10), pady=(5, 10), anchor="nw", fill=tk.BOTH, expand=True)

# Input Row
calendar_input_frame = ttk.Frame(calendar_frame, style='widget.TFrame')
calendar_input_frame.pack(fill=tk.X, padx=5, pady=(5, 10))

min_date_label = ttk.Label(calendar_input_frame, text="Start Date", style='widget.TLabel')
min_date_label.grid(row=0, column=0, sticky="w", padx=(0, 5))
min_date_entry = tk.Entry(calendar_input_frame, width=16)
min_date_entry.grid(row=0, column=1, sticky="w", padx=(0, 15))
add_placeholder(min_date_entry, "YYYY-MM-DD")
min_date_entry.bind("<FocusOut>", on_date_entry_update)
min_date_entry.bind("<Return>", on_date_entry_update)

max_date_label = ttk.Label(calendar_input_frame, text="End Date", style='widget.TLabel')
max_date_label.grid(row=0, column=2, sticky="w", padx=(0, 5))
max_date_entry = tk.Entry(calendar_input_frame, width=16)
max_date_entry.grid(row=0, column=3, sticky="w", padx=(0, 15))
add_placeholder(max_date_entry, "YYYY-MM-DD")
max_date_entry.bind("<FocusOut>", on_date_entry_update)
max_date_entry.bind("<Return>", on_date_entry_update)

calculate_button = ttk.Button(calendar_input_frame, text="Calculate", style ='widget.TButton')
calculate_button.grid(row=0, column=4, sticky="w")
calculate_button.config(command=on_calculate)

# Calendar navigation
calendar_nav_frame = ttk.Frame(calendar_frame, style='widget.TFrame')
calendar_nav_frame.pack(fill=tk.X, padx=5, pady=(0, 5))

prev_month_button = ttk.Button(calendar_nav_frame, text="<", width=3, style='widget.TButton')
prev_month_button.pack(side=tk.LEFT)
prev_month_button.config(command=prev_month)

month_label = ttk.Label(calendar_nav_frame, text="", style='widget_header.TLabel')
month_label.pack(side=tk.LEFT, expand=True)

next_month_button = ttk.Button(calendar_nav_frame, text=">", width=3, style='widget.TButton')
next_month_button.pack(side=tk.LEFT)
next_month_button.config(command=next_month)

today_button = ttk.Button(calendar_nav_frame, text="Today", style='widget.TButton', command=go_to_today)
today_button.pack(side=tk.LEFT, padx=(10, 0))
today_button.config(command=go_to_today)

# Calendar Canvas

calendar_menu = tk.Menu(root, tearoff=0)

canvas = tk.Canvas(calendar_frame, width=7*80, height=7*60, bg="grey25")
canvas.pack(fill=tk.BOTH, expand=True)
canvas.bind("<Configure>", on_canvas_resize)
canvas.bind("<Button-3>", on_calendar_right_click)

now = datetime.now()
current_year = now.year
current_month = now.month

draw_calendar(canvas, current_year, current_month)

root.mainloop()
