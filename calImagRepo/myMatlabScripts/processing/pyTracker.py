# pyTracker.py

def tracker(movie,mouseColor):
    import cv2
    import numpy as np

    # Read video
    video = cv2.VideoCapture(movie)
    frameN = int(video.get(cv2.CAP_PROP_FRAME_COUNT))
    
    position = np.nan * np.zeros([frameN,2]);
    
    # Initialize tracker
    tracker = cv2.createBackgroundSubtractorMOG2()
    for k in range(frameN):
        # Read a new frame
        ok, frame = video.read()
        if not ok:
            break
        
        frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
               
        frameBlur = cv2.GaussianBlur(frame, (21, 21), 0)#Max value foir filter is 21
                
        foreGround = tracker.apply(frameBlur)
                
        tresholdMask = cv2.threshold(frameBlur, int(mouseColor), 255, cv2.THRESH_BINARY)[1] #threshold 80\
               
        foreGround[np.nonzero(tresholdMask)] = 0
                
        Ymask,Xmask=np.nonzero(foreGround)
        
        if 0<len(Xmask) and 0<len(Ymask):
            position[k,0] = np.mean(Xmask)
            position[k,1] = np.mean(Ymask)
                
    return position