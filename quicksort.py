from random import randrange

def quicksort_subarray(nums, start, end):
    if start >= end:
        return 0
    # Partition
    pivot = nums[end]
    partition = start
    for i in range(start, end):
        if nums[i] <= pivot:
            tmp = nums[partition]
            nums[partition] = nums[i]
            nums[i] = tmp
            partition += 1
    tmp = nums[partition]
    nums[partition] = nums[end]
    nums[end] = tmp
    # Recur
    quicksort_subarray(nums, start, partition-1)
    quicksort_subarray(nums, partition+1, end)


def quicksort(nums):
    quicksort_subarray(nums, 0, len(nums)-1)


while True:
    nums = [randrange(100) for _ in range(10)]
    print(nums)
    quicksort(nums)
    print(nums)
    for i in range(len(nums)-1):
        if nums[i] > nums[i+1]:
            print("Not sorted!")
            break
    input()
