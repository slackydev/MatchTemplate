Pure FPC implementation of FFT based template matching.

Template matching is a technique for finding areas of an image that match (are similar) to a template image (smaller image).

The result matrix is normalized to -1..1, and almost identical to that of OpenCV's matchTemplate `CV_TM_CCOEFF_NORMED`