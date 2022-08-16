#ifndef _RANDOM_UTILS_H
#define _RANDOM_UTILS_H

float randf(unsigned int * seed);

#define ISI_PARALLEL
#ifdef ISI_PARALLEL
template <typename T>
T select_random(std::vector<T> const &values, unsigned int *seed) {
  return values[rand_r(seed) % values.size()];
}
#else
template <typename T>
T select_random(std::vector<T> const &values) {
  return values[rand() % values.size()];
}
#endif

template <typename T>
T select_random_determistic(std::vector<T> const &values, std::vector<float> const &weights, float value) {
  if (values.empty()) {
    throw std::invalid_argument("Values list must not be empty.");
  }
  float total = 0.0f;
  for (auto const &w : weights) {
    if (w < 0) {
      throw std::invalid_argument("Weights must not be negative");
    }
    total += w;
  }

  float r = value * total;
  float curr = 0.0f;
  int i = -1;
  while (curr <= r && (i < 0 || i < (int)values.size() - 1)) {
    i++;
    curr += weights[i];
  }
  return values[i];

}

template <typename T>
T select_random(std::vector<T> const &values, std::vector<float> const &weights, unsigned int * seed) {
  return select_random_determistic<T>(values, weights, randf(seed));
}

#endif // _RANDOM_UTILS_H
