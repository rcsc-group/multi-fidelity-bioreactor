"""Trim fully transparent borders of the def_*_light.png panels (the field panels carry letterbox padding)."""
from pathlib import Path

from PIL import Image

for f in sorted(Path(__file__).resolve().parent.parent.joinpath("public").glob("def_*_light.png")):
    im = Image.open(f).convert("RGBA")
    bb = im.getchannel("A").point(lambda a: 255 if a > 8 else 0).getbbox()
    if bb:
        im.crop(bb).save(f)
        print(f.name, im.size, "->", (bb[2] - bb[0], bb[3] - bb[1]))
