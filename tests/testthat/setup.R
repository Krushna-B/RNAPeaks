# Route all graphics-device output to a null device for the duration of the test
# run. Several functions (e.g. draw_plot -> find_left_margin) measure text via
# grid, which opens the default device; headless that default is pdf(), leaving a
# stray Rplots.pdf. Opening pdf(NULL) here absorbs that output and writes nothing.
pdf(NULL)
withr::defer(dev.off(), teardown_env())
