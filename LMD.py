import re

print("Hello World")

expr_object = re.compile("(-?[0-9]*\.[0-9]*)(?:\s+)(-?[0-9]*\.[0-9]*)(?:\s+)(-?[0-9]*\.[0-9]*)(?:\s+)(-?[0-9]*\.[0-9]*)(?:\s+)(-?[0-9]*\.[0-9]*)")
result = expr_object.match("-0000.66      1111.00   -22222.33         3333.11 -4444.6")
print(result.group(0))
print(result.group(1))
print(result.group(2))
print(result.group(3))
print(result.group(4))
# print(result.group(5))
# print(result.group(6))
# print(result.group(7))
