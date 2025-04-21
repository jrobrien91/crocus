import os
from datetime import datetime

import os
import shutil
from datetime import datetime, timedelta

# === CONFIG ===
source_base = "/home/obrienj/crocus/"  # full directory with ALL plots
target_base = "./CROCUS/micronet"      # public directory near index.html
days_back = 4                          # files newer than this get copied
output_file = "micronet_index.html"

# === SETUP ===
now = datetime.now()
cutoff = now - timedelta(days=days_back)

for root, dirs, files in os.walk(source_base):
    for fname in files:
        if not fname.endswith("timeseries.png"):
            continue

        full_path = os.path.join(root, fname)
        mod_time = datetime.fromtimestamp(os.path.getmtime(full_path))

        if mod_time >= cutoff:
            # Calculate relative path from source_base
            rel_path = os.path.relpath(full_path, source_base)
            target_path = os.path.join(target_base, rel_path)

            # Make sure target folder exists
            os.makedirs(os.path.dirname(target_path), exist_ok=True)

            # Copy file
            shutil.copy2(full_path, target_path)
            print(f"Copied: {rel_path}")

# Start building the HTML string
html = """<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Quicklook Plots</title>
    <style>
        body { font-family: Arial, sans-serif; padding: 20px; }
        img.thumb { max-width: 300px; margin: 10px; border: 1px solid #ccc; border-radius: 4px; }
        .plot-block { margin-bottom: 40px; }
        .timestamp { color: #666; font-size: 0.9em; }
        .folder-toggle { cursor: pointer; background: #eee; padding: 8px 12px; border-radius: 5px; margin: 10px 0; }
        .folder-content { display: none; margin-left: 20px; }
        .folder-content.open { display: block; }
    </style>
    <script>
        function toggleFolder(id) {
            const el = document.getElementById(id);
            el.classList.toggle("open");
        }
    </script>
</head>
<body>
<h1>Quicklook Plot Dashboard</h1>
"""

# Build the dynamic content
months = sorted(os.listdir(target_base))
for month in months:
    month_path = os.path.join(target_base, month)
    if not os.path.isdir(month_path):
        continue

    html += f'<div class="plot-block" id="{month}">\n'
    html += f'<div class="folder-toggle" onclick="toggleFolder(\'{month}_content\')">▶ {month}/</div>\n'
    html += f'<div class="folder-content" id="{month}_content">\n'

    for fname in sorted(os.listdir(month_path)):
        if fname.endswith("timeseries.png"):
            fpath = os.path.join(target_base, month, fname)
            timestamp = datetime.fromtimestamp(os.path.getmtime(fpath))
            html += f'''
<div class="plot-item">
    <h4>{fname}</h4>
    <div class="timestamp">Last updated: {timestamp.strftime('%Y-%m-%d %H:%M:%S')}</div>
    <a href="{fpath}" target="_blank"><img class="thumb" src="{fpath}" alt="{fname}"></a>
</div>
'''

    html += "</div></div>\n"

# Close HTML
html += "</body>\n</html>"

# Write to file
with open(output_file, "w") as f:
    f.write(html)

print(f"HTML saved to {output_file}")
