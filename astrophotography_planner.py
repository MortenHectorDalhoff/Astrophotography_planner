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