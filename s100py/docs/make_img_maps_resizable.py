""" sphinx pages with graphviz images for inheritance that have 'maps' to make hotlinks on the image
don't work properly when the image is resized.  This javascript code makes the map be resized to the new dimension of the image on each resize
"""
import glob
import os

from bs4 import BeautifulSoup
from PIL import Image

script = """
    window.onload = function () {{
        var ImageMap = function (orig_width, map, img) {{
                var n,
                    areas = map.getElementsByTagName('area'),
                    len = areas.length,
                    coords = [],
                    previousWidth = orig_width;
                for (n = 0; n < len; n++) {{
                    coords[n] = areas[n].coords.split(',');
                }}
                this.resize = function () {{
                    var n, m, clen,
                        x = img.offsetWidth / previousWidth;
                    for (n = 0; n < len; n++) {{
                        clen = coords[n].length;
                        for (m = 0; m < clen; m++) {{
                            coords[n][m] *= x;
                        }}
                        areas[n].coords = coords[n].join(',');
                    }}
                    previousWidth = img.offsetWidth;
                    return true;
                }};
                window.onresize = this.resize;
            }},
            {image_maps}
        return;
    }}
"""

image_map_template = """
imageMap_{cnt} = new ImageMap({img_width}, document.getElementById('{map_id}'), document.getElementById('{img_id}'));
imageMap_{cnt}.resize();
"""

# iterate all html files
original_dir = os.getcwd()
for fname in glob.glob("./**/*.html", recursive=True):
    print("working on", fname)
    os.chdir(original_dir)
    chg_to_dir, local_fname = os.path.split(fname)
    os.chdir(chg_to_dir)
    # make a list of image templates for this file that need resizng logic
    mappings = []
    soup = BeautifulSoup(open(local_fname, "rb"), 'html.parser')
    # find all <img> elements
    imgs = soup.find_all("img")

    for img_cnt, img in enumerate(imgs):
        # if a map is used to make an image have links
        if "usemap" in img.attrs:
            map_id = img.attrs['usemap'][1:]
            map = soup.find(name="map", id=map_id)
            # give the img an id if it deesn't have one already
            if img.has_attr("id"):
                img_id = img['id']
            else:
                img_id = img.attrs['usemap'][1:] + "_img"
                img.attrs['id'] = img_id

            # get the image size on disk (don't know how to dynamically do this in javascript)
            try:
                im = Image.open(img['src'])
            except Exception as e:
                raise e
            img_width = im.width
            # create a string for this imag/map pair and append it to the list
            mappings.append(image_map_template.format(cnt=img_cnt, img_width=img_width, map_id=map_id, img_id=img_id))
    if mappings:
        # insert the list of img/maps into the script code
        script_to_insert = script.format(image_maps="\n".join(mappings))
        # insert the script into the html
        script_tag = soup.new_tag("script")
        script_tag.string = script_to_insert
        soup.head.insert_before(script_tag)
        with open(local_fname, "w", encoding='utf-8') as output_file:
            output_file.write(str(soup))
            print("found image maps and rewrote", fname)
