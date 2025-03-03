//
// Created by Oscar Holroyd on 04/02/2025.
//

#ifndef PROGRESS_BAR_H
#define PROGRESS_BAR_H

#include <ctime>
#include <string>

enum TimeUnit { SEC, MIN, HR, DAY, YR, KYR, MYR, GYR };

/**
 * @brief Convert a time unit to a string.
 *
 * @param unit The time unit.
 * @return std::string The string representation of the time unit.
 */
std::string time_unit_str(enum TimeUnit unit);

/**
 * @brief Convert a time (with units) to the most appropriate unit.
 *
 * @param time The time.
 * @param unit The time unit.
 */
void time_convert(float *time, enum TimeUnit *unit);

/**
 * @brief Convert a time to a string with a specified format.
 *
 * @param time The time to convert.
 * @param unit The unit of the time.
 * @return std::string The formatted time.
 */
std::string time_format(float time, TimeUnit unit);

class ProgressBar {
public:
  int width;
  int total;
  float progress;

  /**
   * @brief Construct a new Progress Bar object
   *
   * @param total The total number of steps.
   * @param width The width of the progress bar.
   * @param post_print A function that should be called after the progress bar
   * is printed. It should take a single void* argument and print any additional
   * information to the end of the line. Must not print a newline.
   * @param period The minimum time in seconds between updates.
   * @param rank The MPI rank of the process.
   */
  ProgressBar(int total, int width = 50, void (*post_print)(void *) = nullptr,
              float period = 1.0, int rank = 0);

  /**
   * @brief Start the progress bar.
   */
  void start();

  /**
   * @brief Update the progress bar.
   *
   * @param progress The current progress.
   * @param data A pointer to any data that should be passed to the post_print
   * function.
   * @param force If true, force an update even if the period has not elapsed.
   */
  void update(int progress, void *data = nullptr, bool force = false);

  /**
   * @brief End the progress bar.
   */
  void end(void *data = nullptr);

private:
  void (*post_print)(void *);

  float period;
  float last_update;
  std::clock_t start_event;
  std::clock_t mid_event;
  float elapsed_time;
  int rank;

  /**
   * @brief Update the progress bar.
   */
  void update_bar();
};

#endif // PROGRESS_BAR_H
