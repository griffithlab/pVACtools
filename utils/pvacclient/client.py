import os
import site
import http.server
import socketserver
#import multiprocessing
#import webbrowser
#from subprocess import Popen

def main():
	client_dir = ''
	site_dirs = site.getsitepackages()
	for path in site_dirs:
		tmp_path = os.path.join(
			path,
			'pvacseq-client'
		)
		if os.path.isdir(tmp_path):
			client_dir = tmp_path
			break
	if not len('client_dir'):
		sys.exit("Unable to locate the frontend!")

	print("Launching Frontend Server")
	os.chdir(client_dir)
	Handler = http.server.SimpleHTTPRequestHandler
	httpd = socketserver.TCPServer(("localhost",8000), Handler)
	# Simplest way to just serve the client on http://localhost:8000. Call `pvacseq-api --nogui` in another window
	# if you want it to connect to the API, or just `pvacseq-api` if the code that calls up the client ever gets removed
	# from the initialization of the API (see last comment section). Exiting out can get a bit messy though.
	httpd.serve_forever()

	# If you want to automatically open localhost:8000 in a new tab, use the command below.
	# You won't have to Ctrl+C to close the server since it doesn't take control of what
	# the console prints out like serve_forever() does.
	#
	# webbrowser.open_new_tab('http://localhost:8000')

	# If you want to start the API after starting the front end without calling pvacseq-api in a separate window,
	# use this instead of httpd.serve_forever() or webbrowser.open_new_tab() so that the code
	# will continue executing even after starting the client. You might not need to store it in the server_process
	# variable and be able to just call multiprocessing.Process()
	#
	# server_process = multiprocessing.Process(target=httpd.serve_forever).start()

	# One method of starting the API, but it won't go through the cleanup process it would if you just called
	# pvacseq-api and Ctrl-C'd to exit, so make sure to manually kill processes on ports 8000, 8080, and 5006.
	# Make sure to remove the section starting with "if '--nogui' not in args:" in the initialize() function
	# of pvacapi/controllers/utils.py if you don't want the API to start the client.
	# You might not need to store this in the api variable and just be able to call Popen()
	#
	# api = Popen(
	# 	[
	# 		sys.executable,
	# 		os.path.join(os.path.dirname(os.path.dirname(__file__)),'pvacapi', 'app.py'),
	# 		'api'
	# 	]
	# )

if __name__ == '__main__':
	main()