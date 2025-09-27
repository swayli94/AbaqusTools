# -*- coding: utf-8 -*-
"""
Abaqus/CAE 2023: create a 'plasma' Spectrum.
"""

# -------- User settings --------
SPECTRUM_NAME = "plasma"
N_COLORS = 64           # 64-256 typical
# --------------------------------

from abaqus import session

# ---------- Color helpers ----------
def rgb01_to_hex(rgb):
    r = int(round(max(0.0, min(1.0, float(rgb[0]))) * 255))
    g = int(round(max(0.0, min(1.0, float(rgb[1]))) * 255))
    b = int(round(max(0.0, min(1.0, float(rgb[2]))) * 255))
    return "#{:02X}{:02X}{:02X}".format(r, g, b)

def get_plasma_hex_colors(n=256):
    try:
        import matplotlib.cm as cm
        cmap = cm.get_cmap("plasma", n)
        return [rgb01_to_hex(cmap(i)[:3]) for i in range(cmap.N)]
    except Exception:
        # Fallback key plasma colors (approximate), with linear interpolation
        key_hex = [
            "#0D0887","#3A049A","#5C01A6","#7E03A8","#9C179E",
            "#B52F8C","#CC4778","#DD5E66","#EA7851","#F3943E",
            "#FBB227","#F7D13D","#F0F921"
        ]
        def hex_to_rgb01(h):
            h = h.lstrip("#")
            return (int(h[0:2],16)/255.0, int(h[2:4],16)/255.0, int(h[4:6],16)/255.0)
        keys = [hex_to_rgb01(h) for h in key_hex]

        import math
        out = []
        m = len(keys) - 1
        for i in range(n):
            t = float(i) / (n - 1 if n > 1 else 1.0)
            pos = t * m
            a = int(math.floor(pos))
            b = min(a + 1, m)
            lt = pos - a
            r = keys[a][0]*(1-lt) + keys[b][0]*lt
            g = keys[a][1]*(1-lt) + keys[b][1]*lt
            b_ = keys[a][2]*(1-lt) + keys[b][2]*lt
            out.append(rgb01_to_hex((r, g, b_)))
        return out

# ---------- Spectrum management ----------
def delete_spectrum_if_exists(name):
    if name in session.spectrums.keys():
        try:
            del session.spectrums[name]
        except:
            pass

def create_spectrum(name, color_hex_list):
    colors_tuple = tuple(str(c) for c in color_hex_list)
    spec = session.Spectrum(name=name, colors=colors_tuple)
    return spec

# ---------- Main ----------
def main():
    colors_hex = get_plasma_hex_colors(N_COLORS)
    delete_spectrum_if_exists(SPECTRUM_NAME)
    spec = create_spectrum(SPECTRUM_NAME, colors_hex)
    print("Spectrum '{}' created with {} colors.".format(spec.name, len(colors_hex)))

if __name__ == "__main__":
    main()
