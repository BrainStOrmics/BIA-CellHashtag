# Message queues

[Frontend](/oss/python/langchain/frontend/overview)[Patterns](/oss/python/langchain/frontend/markdown-messages)

# Message queues
Copy page

Queue multiple messages and manage them while the agent processes sequentiallyCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Message queuing lets users send multiple messages in rapid succession without waiting for the agent to finish processing the current one. Each message is enqueued server-side and processed sequentially, giving you full visibility and control over the pending queue.

This feature requires the [LangGraph Agent Server](/langsmith/local-server). Run your agent locally with `langgraph dev` or [deploy it to LangSmith](/langsmith/deployment) to use this pattern.

## [​](#why-message-queues)Why message queues?

In a typical chat interface, users must wait for the agent to finish responding before sending another message. This creates friction in several scenarios:

{' ' * (self.list_depth - 1)}- **Batch questions**: a user wants to ask five related questions at once rather than waiting for each answer

{' ' * (self.list_depth - 1)}- **Follow-up chains**: submitting clarifications or additional context while the agent is still working

{' ' * (self.list_depth - 1)}- **Automated testing sequences**: programmatically sending a series of prompts to validate agent behavior

{' ' * (self.list_depth - 1)}- **Data entry workflows**: feeding structured inputs one after another for processing

Message queuing solves this by accepting all submissions immediately and processing them in order.

## [​](#how-it-works)How it works

Under the hood, LangGraph uses `multitaskStrategy: "enqueue"` to manage concurrent submissions. When a message is submitted while the agent is already processing, it gets added to a server-side queue. Once the current run completes, the next queued message is picked up automatically.
The `useStream` hook exposes a `queue` property that provides real-time visibility into pending messages:

| 
 | Property | Type | Description
 | `queue.entries` | `QueueEntry[]` | Array of all pending queue entries
 | `queue.size` | `number` | Number of entries currently in the queue
 | `queue.cancel(id)` | `(id: string) => Promise<void>` | Cancel a specific queued entry by ID
 | `queue.clear()` | `() => Promise<void>` | Cancel all queued entries
Each `QueueEntry` object contains:

| 
 | Field | Type | Description
 | `id` | `string` | Unique identifier for this queue entry
 | `values` | `object` | The input values (including messages) that were submitted
 | `options` | `object` | Any additional options passed with the submission
 | `createdAt` | `string` | ISO timestamp of when the entry was created

## [​](#setting-up-usestream)Setting up `useStream`

Define a TypeScript interface matching your agent’s state schema and pass it as a type parameter to `useStream` for type-safe access to state values. In the examples below, replace `typeof myAgent` with your interface name:

```
import type { BaseMessage } from "@langchain/core/messages";

interface AgentState {
 messages: BaseMessage[];
}

```

ReactVueSvelteAngular
```
import { useStream } from "@langchain/react";

function Chat() {
 const stream = useStream<typeof myAgent>({
 apiUrl: "http://localhost:2024",
 assistantId: "message_queue",
 });

 const handleSubmit = (text: string) => {
 stream.submit({
 messages: [{ type: "human", content: text }],
 });
 };

 // Access queue state
 const pendingCount = stream.queue.size;
 const entries = stream.queue.entries;

 return (
 <div>
 <MessageList messages={stream.messages} />
 {pendingCount > 0 && <QueueList entries={entries} queue={stream.queue} />}
 <ChatInput onSubmit={handleSubmit} />
 </div>
 );
}

```

## [​](#displaying-the-queue)Displaying the queue

Build a `QueueList` component that shows each pending message with a cancel button. This gives users visibility into what’s waiting and the ability to remove items they no longer need.

```
function QueueList({ entries, queue }) {
 return (
 <div className="queue-panel">
 <div className="queue-header">
 <span>Queued messages ({entries.length})</span>
 <button onClick={() => queue.clear()}>Clear all</button>
 </div>
 <ul className="queue-entries">
 {entries.map((entry) => {
 const text = entry.values?.messages?.[0]?.content ?? "Unknown";
 return (
 <li key={entry.id} className="queue-entry">
 <span className="queue-text">{text}</span>
 <span className="queue-time">
 {new Date(entry.createdAt).toLocaleTimeString()}
 </span>
 <button
 className="queue-cancel"
 onClick={() => queue.cancel(entry.id)}
 >
 Cancel
 </button>
 </li>
 );
 })}
 </ul>
 </div>
 );
}

```

Display the first few characters of each queued message as a preview so users can quickly identify which items to cancel without reading full messages.

## [​](#cancelling-queued-messages)Cancelling queued messages

You have two levels of cancellation:

### [​](#cancel-a-single-entry)Cancel a single entry

Remove a specific message from the queue by its ID. The agent will skip it and move to the next entry.

```
await queue.cancel(entryId);

```

### [​](#clear-the-entire-queue)Clear the entire queue

Remove all pending messages at once. Useful when the user changes context or wants to start over.

```
await queue.clear();

```

Cancelling a queue entry only affects messages that have **not yet started
processing**. If the agent is already working on a message, cancelling it from
the queue has no effect. Use `stream.stop()` to interrupt the current run.

## [​](#chaining-follow-up-submissions-with-oncreated)Chaining follow-up submissions with `onCreated`

The `onCreated` callback fires when a new run is created, giving you a hook to submit follow-up messages programmatically. This is useful for building multi-step workflows where the next question depends on the previous submission being accepted.

```
stream.submit(
 { messages: [{ type: "human", content: "What is quantum computing?" }] },
 {
 onCreated(run) {
 console.log("Run created:", run.run_id);
 // Chain a follow-up
 stream.submit({
 messages: [{ type: "human", content: "Give me a simple analogy." }],
 });
 },
 }
);

```

This pattern naturally fills the queue. The first message starts processing
immediately, and the follow-up is queued behind it.

## [​](#starting-a-new-thread)Starting a new thread

When a user wants to begin a fresh conversation, use `switchThread(null)` to create a new thread. This clears the current message history and queue.
ReactVueSvelteAngular
```
function NewThreadButton() {
 const stream = useStream<typeof myAgent>({ /* ... */ });

 return (
 <button onClick={() => stream.switchThread(null)}>
 New conversation
 </button>
 );
}

```

## [​](#complete-example)Complete example

Putting it all together, here is a full chat component with queue management:

```
function QueueChat() {
 const stream = useStream<typeof myAgent>({
 apiUrl: "http://localhost:2024",
 assistantId: "message_queue",
 });

 const [input, setInput] = useState("");

 const handleSubmit = () => {
 if (!input.trim()) return;
 stream.submit({
 messages: [{ type: "human", content: input.trim() }],
 });
 setInput("");
 };

 return (
 <div className="chat-container">
 <header>
 <h2>Queue Chat</h2>
 <button onClick={() => stream.switchThread(null)}>New thread</button>
 </header>

 <div className="messages">
 {stream.messages.map((msg, i) => (
 <MessageBubble key={i} message={msg} />
 ))}
 {stream.isLoading && <TypingIndicator />}
 </div>

 {stream.queue.size > 0 && (
 <div className="queue-panel">
 <strong>Queued ({stream.queue.size})</strong>
 <button onClick={() => stream.queue.clear()}>Clear all</button>
 {stream.queue.entries.map((entry) => (
 <div key={entry.id} className="queue-item">
 <span>{entry.values?.messages?.[0]?.content}</span>
 <button onClick={() => stream.queue.cancel(entry.id)}>×</button>
 </div>
 ))}
 </div>
 )}

 <form onSubmit={(e) => { e.preventDefault(); handleSubmit(); }}>
 <input
 value={input}
 onChange={(e) => setInput(e.target.value)}
 placeholder="Type a message (you can send multiple!)"
 />
 <button type="submit">Send</button>
 </form>
 </div>
 );
}

```

## [​](#best-practices)Best practices

{' ' * (self.list_depth - 1)}- **Limit queue size**: While there is no hard client-side limit on queue size,
be mindful that very large queues can degrade user experience. Consider
showing a warning when the queue exceeds a reasonable threshold (e.g., 10
items).

{' ' * (self.list_depth - 1)}- **Show queue position**: Number each queued item so users know the processing order.

{' ' * (self.list_depth - 1)}- **Preserve input focus**: Keep the input field focused after submission so users can type the next message immediately.

{' ' * (self.list_depth - 1)}- **Animate transitions**: Smoothly move items from the queue panel into the message list as they start processing.

{' ' * (self.list_depth - 1)}- **Handle errors gracefully**: If a queued message fails, surface the error without blocking subsequent queue entries.

{' ' * (self.list_depth - 1)}- **Debounce rapid submissions**: For automated or programmatic submissions, add a small delay between messages to avoid overwhelming the server.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/frontend/message-queues.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Structured outputPrevious](/oss/python/langchain/frontend/structured-output)[Join & rejoin streamsNext](/oss/python/langchain/frontend/join-rejoin)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
