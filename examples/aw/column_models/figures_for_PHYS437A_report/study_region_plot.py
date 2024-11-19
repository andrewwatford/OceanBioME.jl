import matplotlib.pyplot as plt
import cartopy.crs as ccrs


def main():
    lat_min, lat_max = 45, 50
    lon_min, lon_max = -35, -25
    x, y = [lon_min, lon_min, lon_max, lon_max, lon_min], [lat_min, lat_max, lat_max, lat_min, lat_min]
    extent = [-80, 20, 0, 80]

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    # make the map global rather than have it zoom in to
    # the extents of any plotted data
    ax.set_global()

    ax.stock_img()
    ax.coastlines()
    ax.gridlines(draw_labels=True)

    ax.plot(x, y, c = 'k', linewidth=2, transform=ccrs.Geodetic())
    ax.fill(list(reversed(x)), list(reversed(y)), color = 'black', alpha = 0.4, transform=ccrs.Geodetic())
    ax.set_extent(extent, crs=ccrs.Geodetic())
    plt.savefig("./study_region.png", dpi=400)


    plt.show()


if __name__ == '__main__':
    main()