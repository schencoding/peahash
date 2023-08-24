//
// This file is to generate unique keys from random keys
//
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <unordered_set>
#include <iostream>
#include <fstream>

#define KEY_LEN 8
#define SAMPLE_NUM  5000000L
#define MAX_LEN     250000000L

std::string insert_pth = "../keys/com-orkut.ungraph";

std::string orig_pth = insert_pth + ".txt";
std::string skew_sp_pth = insert_pth + "_sp5M";
std::string uniq_sp_pth = insert_pth + "_uniq_sp";

void shuffle(uint64_t *workload, uint64_t number) {
  srand(time(nullptr));
  for (size_t i = 0; i < 2*number; i++)
  {
    int replace_id = rand() % number;
    uint64_t swap_one = workload[replace_id];
    workload[replace_id] = workload[0];
    workload[0] = swap_one;
  }
}

void resample(const uint64_t *workload, uint64_t *sample_arr, uint64_t total_num, uint64_t sample_num) {
  srand(time(nullptr));
  for (size_t i = 0; i < sample_num; i++) {
    int resample_id = rand() % total_num;
    sample_arr[i] = workload[resample_id];
  }
}

// before: malloc now
// return: the distinct elements' num
size_t distinctSample(uint64_t *orig, uint64_t *now, size_t cap, size_t sample_num){
  std::unordered_set<uint64_t> table;
  size_t j = 0;
  for (size_t i = 0; i < cap; ++i) {
    if (table.count(orig[i]) == 0){
      table.emplace(orig[i]);
      now[j] = orig[i];
      j++;
    }
    if (j == sample_num){
      break;
    }
  }
  return j;
}


int main() {
  auto *array = static_cast<uint64_t *>(malloc(MAX_LEN * KEY_LEN));

  std::ifstream infile;
  infile.open(orig_pth.data(), std::ios::in);
  if (!infile.is_open()) {
      std::cout << "读取文件失败" << std::endl;
      return 0;
  }
  char buf[1024] = {0};
  uint64_t num;
  uint64_t idx = 0;
  while (infile >> buf) {
    num = std::strtoull(buf, (char**) nullptr, 10);
    array[idx++] = num;
  }
  infile.close();
  uint64_t insert_num = idx;
  shuffle(array, insert_num);
  FILE *fp = fopen(insert_pth.data(),"w");
  fwrite(array, KEY_LEN, insert_num, fp);
  fclose(fp);

  auto *sample_array = static_cast<uint64_t *>(malloc(SAMPLE_NUM * KEY_LEN));
  resample(array, sample_array, insert_num, SAMPLE_NUM);
  fp = fopen(skew_sp_pth.data(),"w");
  fwrite(sample_array, KEY_LEN, SAMPLE_NUM, fp);
  fclose(fp);

  shuffle(array, insert_num);
  size_t sample_size = distinctSample(array, sample_array, insert_num, SAMPLE_NUM);
  std::cout << "The sample num: " << sample_size << std::endl;
  fp = fopen(uniq_sp_pth.data(),"w");
  fwrite(sample_array, KEY_LEN, sample_size, fp);
  fclose(fp);

  free(array);
  free(sample_array);

  return 0;
}

