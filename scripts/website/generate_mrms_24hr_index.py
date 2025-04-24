import os
from datetime import datetime

import os
import shutil
from datetime import datetime, timedelta

# === CONFIG ===
source_base = "/home/obrienj/crocus-analysis/mrms/qpe-24/"  # full directory with ALL plots
target_base = "/home/obrienj/public_html/CROCUS/mrms-24hr-qpe/"      # public directory near index.html
html_target = "./CROCUS/mrms-24hr-qpe/"
days_back = 3                          # files newer than this get copied
remove_days_back = 4
output_file = "/home/obrienj/public_html/mrms_qpe_index.html"

# === SETUP ===
now = datetime.now()
cutoff = now - timedelta(days=days_back)
remove_cutoff = now - timedelta(days=remove_days_back)

for root, dirs, files in os.walk(source_base):
    for fname in files:
        if not fname.startswith("mrms-24hr-qpe"):
            print('hey startswith works')
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

# === DELETE OLD FILES FROM TARGET ===
for root, dirs, files in os.walk(target_base):
    for fname in files:
        if not fname.startswith("mrms-24hr-qpe"):
            continue

        target_path = os.path.join(root, fname)
        mod_time = datetime.fromtimestamp(os.path.getmtime(target_path))

        if mod_time < remove_cutoff:
            os.remove(target_path)
            rel_path = os.path.relpath(target_path, target_base)

            print(f"Deleted (too old): {rel_path}")


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
<h1>MRMS 24hr Quantitative Precipitation Estimate (QPE) with CROCUS Micronet Observations</h1>
"""

# Build the dynamic content
##html += f'<div class="plot-block" id="CROCUS">\n'
##html += f'<div class="folder-toggle" onclick="toggleFolder(CROCUS)">▶ "MRMS 24hr QPE - CROCUS"/</div>\n'
##html += f'<div class="folder-content" id="MRMS 24hr QPE - CROCUS">\n'

for fname in sorted(os.listdir(target_base)):
    if fname.startswith("mrms-24hr-qpe"):
        fpath = os.path.join(target_base, fname)
        timestamp = datetime.fromtimestamp(os.path.getmtime(fpath))
        hpath = os.path.join(html_target, fname)
        html += f'''
<div class="plot-item">
    <h4>{fname}</h4>
    <div class="timestamp">Last updated: {timestamp.strftime('%Y-%m-%d %H:%M:%S')}</div>
    <a href="{hpath}" target="_blank"><img class="thumb" src="{hpath}" alt="{fname}"></a>
</div>
'''

html += "</div></div>\n"

# Close HTML
html += "</body>\n</html>"

# Write to file
with open(output_file, "w") as f:
    f.write(html)

print(f"HTML saved to {output_file}")
