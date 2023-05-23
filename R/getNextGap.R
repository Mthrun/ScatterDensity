getNextGap = function(data, currentIndex, minIndex, maxIndex, direction) {
  nextGap = 0
  pointer = currentIndex
  while (pointer > minIndex && pointer < maxIndex) { # TURE
    pointer = pointer + direction
    nextGap = abs(data[pointer] - data[currentIndex])
    if (nextGap > 0) { # Ignores duplicates in data
      return(pointer)
    }
  }
  return(pointer)
}
