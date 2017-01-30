import math

# Returns a percentile of an array. It is assumed that the array is already
# sorted.
#
# The array must not be empty.
#
# Percentile must have a value between [0, 1.0].
def percentile(array, percentile):
  assert len(array)
  assert percentile >= 0.0 and percentile <= 1.0
  index = int(math.ceil(percentile * len(array))) - 1
  return array[index]
