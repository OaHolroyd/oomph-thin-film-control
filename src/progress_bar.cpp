//
// Created by Oscar Holroyd on 04/02/2025.
//

#include "progress_bar.h"

#include <cstdio>
#include <iomanip>
#include <cstdarg>


std::string time_unit_str(TimeUnit unit) {
  switch (unit) {
    case SEC:
      return "s";
    case MIN:
      return "min";
    case HR:
      return "hr";
    case DAY:
      return "day";
    case YR:
      return "yr";
    case KYR:
      return "kyr";
    case MYR:
      return "Myr";
    case GYR:
      return "Gyr";
  }
  return "";
}

/**
 * @brief Convert a time to a larger unit.
 *
 * @param time The time to convert.
 * @param unit The unit of the time.
 */
void time_convert_up(float *time, TimeUnit *unit) {
  switch (*unit) {
    case SEC:
      *time /= 60.0;
      *unit = MIN;
      break;
    case MIN:
      *time /= 60.0;
      *unit = HR;
      break;
    case HR:
      *time /= 24.0;
      *unit = DAY;
      break;
    case DAY:
      *time /= 365.25;
      *unit = YR;
      break;
    case YR:
      *time /= 1000.0;
      *unit = KYR;
      break;
    case KYR:
      *time /= 1000.0;
      *unit = MYR;
      break;
    case MYR:
      *time /= 1000.0;
      *unit = GYR;
      break;
    default:
      // do nothing
      break;
  }
}


/**
 * @brief Convert a time to a smaller unit.
 *
 * @param time The time to convert.
 * @param unit The unit of the time.
 */
void time_convert_down(float *time, TimeUnit *unit) {
  switch (*unit) {
    case MIN:
      *time *= 60.0;
      *unit = SEC;
      break;
    case HR:
      *time *= 60.0;
      *unit = MIN;
      break;
    case DAY:
      *time *= 24.0;
      *unit = HR;
      break;
    case YR:
      *time *= 365.25;
      *unit = DAY;
      break;
    case KYR:
      *time *= 1000.0;
      *unit = YR;
      break;
    case MYR:
      *time *= 1000.0;
      *unit = KYR;
      break;
    case GYR:
      *time *= 1000.0;
      *unit = MYR;
      break;
    default:
      // do nothing
      break;
  }
}


void time_convert(float *time, TimeUnit *unit) {
  // change to larger units until time is less than 1.0
  while (*time > 1.0 && *unit != GYR) {
    time_convert_up(time, unit);
  }

  // change back down to the first unit where time that is greater than 1.0
  while (*time < 1.0 && *unit != SEC) {
    time_convert_down(time, unit);
  }
}


std::string format_string(const char *fmt, ...) {
  const int len = 512;
  std::string s(len, '\0');
  va_list args;
  va_start(args, fmt);
  auto r = std::vsnprintf(&s[0], len, fmt, args);
  va_end(args);
  if (r < 0) return "??";
  s.resize(std::min(r, len));
  return s;
}


std::string time_format(float time, TimeUnit unit) {
  /* handle special cases */
  if (time == 0.0) {
    return "0" + time_unit_str(unit);
  }
  if (time < 0.0) {
    return "-" + time_format(-time, unit);
  }

  /* convert time to the appropriate unit */
  time_convert(&time, &unit);

  /* format the time */
  return format_string("%.1f", time) + time_unit_str(unit);
}


ProgressBar::ProgressBar(int total, int width, void (*post_print)(void *), float period) {
  this->total = total;
  this->progress = 0.0;
  this->period = period;

  this->start_event = 0;
  this->mid_event = 0;

  this->width = width;
  this->post_print = post_print;
}

void ProgressBar::start() {
  this->start_event = std::clock();
  this->progress = 0.0;
  this->last_update = 0.0;
}

void ProgressBar::update(int progress, void *data, bool force) {
  // calculate progress
  this->progress = static_cast<float>(progress) / static_cast<float>(this->total);

  // compute elapsed time
  this->mid_event = std::clock();
  this->elapsed_time = 1000.0 * (float)(this->mid_event - this->start_event) / (float) CLOCKS_PER_SEC;

  // decide whether to print
  if ((this->elapsed_time - this->last_update < this->period * 1000.0) && !force) {
    return;
  }
  this->last_update = this->elapsed_time;

  this->update_bar();
  if (this->post_print != nullptr) {
    this->post_print(data);
  }

  // flush
  fflush(stderr);
}

void ProgressBar::end(void *data) {
  // print update
  this->update(this->total, data, true);

  // print end info
  fprintf(stderr, "\n");
  fprintf(stderr, "Finished. Elapsed time: %3.2fs\n", this->elapsed_time / 1000.0);
}

void ProgressBar::update_bar() {
  int reps = static_cast<int>(this->width * this->progress);

  // restart line
  fprintf(stderr, "\r");

  // percentage complete
  fprintf(stderr, "%3.0f%%", this->progress * 100.0);

  // progress bar
  fprintf(stderr, " [");
  for (int i = 0; i < reps; i++) {
    fprintf(stderr, "|");
  }
  for (int i = reps; i < this->width; i++) {
    fprintf(stderr, " ");
  }
  fprintf(stderr, "]");

  // elapsed time
  fprintf(stderr, " elap: %s", time_format(this->elapsed_time / 1000.0f, SEC).c_str());

  // estimated time remaining
  float total_time = this->elapsed_time / this->progress;
  float remaining_time = total_time - this->elapsed_time;
  fprintf(stderr, ", rem: %s", time_format(remaining_time / 1000.0f, SEC).c_str());
}