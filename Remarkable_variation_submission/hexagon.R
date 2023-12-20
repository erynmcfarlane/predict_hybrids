### hexagon function

hexagon <- function(x_center=NA, y_center=NA, scale=NA, plot_xmin=NA, plot_xmax=NA, plot_ymin=NA, plot_ymax=NA, color=NA)
  ## function allows you to plot fillable hexagon points in plots
  ## x_center and y_center mark the center of the point
  ## scale controls the size of the point relative to the size of the plotting area (e.g., 0.1)
  ## plot_xmin, plot_xmax, plot_ymin, and plot_ymax need to be the xlims and ylims of your plot
  ## color is any color you'd like to use
{
  x_scale <- scale * (plot_xmax - plot_xmin)
  y_scale <- scale * (plot_ymax - plot_ymin)
  ax <- x_center + x_scale
  ay <- y_center
  bx <- x_center + x_scale / 2
  by <- y_center + sqrt(3) * y_scale / 2
  cx <- x_center - x_scale / 2
  cy <- y_center + sqrt(3) * y_scale / 2
  dx <- x_center - x_scale
  dy <- y_center
  ex <- x_center - x_scale / 2
  ey <- y_center - sqrt(3) * y_scale / 2
  fx <- x_center + x_scale / 2
  fy <- y_center - sqrt(3) * y_scale / 2
  x_points <- c(ax, bx, cx, dx, ex, fx)
  y_points <- c(ay, by, cy, dy, ey, fy)
  polygon(x_points, y_points, col=color, border="black")
}
