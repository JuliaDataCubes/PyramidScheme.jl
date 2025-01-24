module PyramidSchemeTylerExt
import Tyler
import PyramidScheme: PyramidProvider, selectlevel_tile, fetch_tile
import MapTiles: Tile, extent
import Colors: RGB
import ColorSchemes: colorschemes
import DimensionalData: DimArray, Near
import TileProviders

Tyler.fetch_tile(p::PyramidProvider, tile::Tile) = fetch_tile(p,tile)
TileProviders.max_zoom(p::PyramidProvider) = p.max_zoom
TileProviders.min_zoom(p::PyramidProvider) = p.min_zoom

end