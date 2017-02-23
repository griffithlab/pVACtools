import watchdog.events
import watchdog.observers
import threading

class pvacseqHandler(watchdog.events.FileSystemEventHandler):

    def __init__(self, *args, **kwargs):
        super().__init__()
        self.subscribers = []

    def on_any_event(self, event):
        threading.Thread(target=pvacseqHandler._runner, args=(self, event), daemon=True).start()

    def subscribe(self, fn, eventType = None):
        self.subscribers.append((eventType, fn))
        return lambda:self.subscribers.remove((eventType, fn)) #return unsubscribe function

    def _runner(self, event):
        for subscriber in self.subscribers:
            if (not subscriber[0]) or isinstance(event, subscriber[0]):
                subscriber[1](event)

def Observe(path):
    observer = watchdog.observers.Observer()
    handler = pvacseqHandler()
    observer.schedule(handler, path, True)
    observer.subscribe = handler.subscribe
    observer.start()
    return observer
