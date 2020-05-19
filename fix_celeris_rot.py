import pyautogui

pyautogui.PAUSE = 1
pyautogui.FAILSAFE = True

size = pyautogui.size()
print(size)
w, h = pyautogui.size()

drag_distance = 1

pyautogui.moveTo(w/2+300, h/2, duration=.25)
pyautogui.dragRel(drag_distance, 0, button='right', duration=.25)